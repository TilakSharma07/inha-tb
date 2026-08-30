"""Step 8 - dock the prepared ligand library into InhA (PDB 4TZK) with AutoDock Vina.

Design notes, both learned the hard way on this machine:

* **Resume is mandatory.** A ~500-ligand screen runs for hours on CPU, and this one was
  killed twice mid-run (once when the launching cell died, once under memory pressure while
  QSAR training competed for the same 8 cores). Vina writes its score into the pose file's
  `REMARK VINA RESULT` line, so a completed ligand is fully recoverable from disk: this
  script reads back every existing pose file and only docks what is missing.
* **Results are written incrementally**, one JSON line per ligand, as each finishes. The
  earlier version accumulated results in memory and wrote once at the end, so a crash at 45%
  discarded 45% of the compute even though the pose files had survived.
* **Worker count times cpu-per-worker must not exceed the core count.** Each worker holds
  its own copy of the affinity grid maps, so oversubscribing costs RAM as well as context
  switches.
"""
import glob, json, os, re, time
import pandas as pd
import multiprocessing as mp
import numpy as np

# script lives in src/, ligands and receptor live in dock/
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dock")
NPROC = max(1, min(4, (os.cpu_count() or 4) // 2))
VINA_CPU, EXH, N_POSES = 2, 12, 9
SCORES_JSONL = f"{D}/docking_scores.jsonl"
RE_SCORE = re.compile(r"REMARK VINA RESULT:\s*(-?\d+\.\d+)")


def prep_one(sdf):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    lid = os.path.basename(sdf)[:-4]
    out = f"{D}/pdbqt/{lid}.pdbqt"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return lid
    m = Chem.MolFromMolFile(sdf, removeHs=False)
    if m is None:
        return None
    try:
        setups = MoleculePreparation().prepare(m)
        txt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None
        open(out, "w").write(txt)
        return lid
    except Exception:
        return None


def scores_from_pose(lid):
    """Recover a finished ligand's energies from its pose file (crash-resume path)."""
    p = f"{D}/poses/{lid}.pdbqt"
    if not (os.path.exists(p) and os.path.getsize(p) > 0):
        return None
    e = [float(x) for x in RE_SCORE.findall(open(p).read())]
    if not e:
        return None
    return dict(lig_id=lid, vina_kcal=e[0], vina_mean_top3=float(np.mean(e[:3])),
                pose_spread=float(np.std(e)), n_poses=len(e), source="pose_file")


_V = None


def _init_worker():
    """One Vina instance per worker: the grid maps are computed once, not per ligand."""
    global _V
    from vina import Vina
    box = json.load(open(f"{D}/box.json"))
    _V = Vina(sf_name="vina", cpu=VINA_CPU, seed=42, verbosity=0)
    _V.set_receptor(f"{D}/receptor.pdbqt")
    _V.compute_vina_maps(center=box["center"], box_size=box["size"])


def dock_one(lid):
    try:
        _V.set_ligand_from_file(f"{D}/pdbqt/{lid}.pdbqt")
        _V.dock(exhaustiveness=EXH, n_poses=N_POSES)
        en = _V.energies(n_poses=N_POSES)
        _V.write_poses(f"{D}/poses/{lid}.pdbqt", n_poses=3, overwrite=True)
        e = [float(x[0]) for x in en]
        return dict(lig_id=lid, vina_kcal=e[0], vina_mean_top3=float(np.mean(e[:3])),
                    pose_spread=float(np.std(e)), n_poses=len(e), source="docked")
    except Exception as exc:
        return dict(lig_id=lid, vina_kcal=None, error=type(exc).__name__)


def main():
    os.makedirs(f"{D}/pdbqt", exist_ok=True)
    os.makedirs(f"{D}/poses", exist_ok=True)
    # dock/library.csv is the selection of record. Globbing ligands/ instead would dock
    # conformer files left behind by an earlier selection - work that is thrown away at
    # the join, because the analysis only ever reads molecules present in the library.
    lib = pd.read_csv(f"{D}/library.csv")
    sdfs = [f"{D}/ligands/{i}.sdf" for i in lib.lig_id]
    missing = [p for p in sdfs if not os.path.exists(p)]
    if missing:
        raise SystemExit(f"{len(missing)} library molecules have no conformer "
                         f"(run src/06_build_library.py first): {missing[:3]}")
    t0 = time.time()
    with mp.Pool(NPROC) as pool:
        ids = [x for x in pool.map(prep_one, sdfs) if x]
    print(f"prepared {len(ids)}/{len(sdfs)} pdbqt in {time.time()-t0:.0f}s", flush=True)

    done, todo = {}, []
    for lid in ids:
        r = scores_from_pose(lid)
        (done.setdefault(lid, r) if r else todo.append(lid))
    print(f"resume: {len(done)} already docked, {len(todo)} to run, "
          f"{NPROC} workers x {VINA_CPU} cpu", flush=True)

    # rewrite the ledger from what is on disk, then append as ligands finish
    with open(SCORES_JSONL, "w") as fh:
        for r in done.values():
            fh.write(json.dumps(r) + "\n")

    t1, n = time.time(), 0
    if todo:
        with mp.Pool(NPROC, initializer=_init_worker) as pool, open(SCORES_JSONL, "a") as fh:
            for r in pool.imap_unordered(dock_one, todo, chunksize=1):
                fh.write(json.dumps(r) + "\n")
                fh.flush()
                n += 1
                if n % 10 == 0:
                    rate = n / max(1e-9, (time.time() - t1) / 60)
                    print(f"  {n}/{len(todo)} docked ({rate:.1f}/min, "
                          f"ETA {(len(todo)-n)/max(rate,1e-9):.0f} min)", flush=True)

    rows = [json.loads(l) for l in open(SCORES_JSONL)]
    json.dump(rows, open(f"{D}/docking_scores.json", "w"))
    good = [r["vina_kcal"] for r in rows if r.get("vina_kcal") is not None]
    print(f"DONE {len(good)}/{len(ids)} scored "
          f"({n} this run, {(time.time()-t1)/60:.1f} min)", flush=True)
    print(f"vina: min={min(good):.2f} median={np.median(good):.2f} "
          f"max={max(good):.2f} kcal/mol", flush=True)


if __name__ == "__main__":
    main()


def prune_orphans(dry_run=True):
    """Report (or delete) conformer/pose files for molecules no longer in the library.

    Re-selecting the library leaves the previous selection's files on disk. They are inert -
    every downstream step joins on dock/library.csv, so an orphan pose is never read - but
    they make `ls dock/poses | wc -l` disagree with the reported screen size, which is the
    kind of discrepancy that costs an hour to re-derive later. Default is dry-run.
    """
    lib = set(pd.read_csv(f"{D}/library.csv").lig_id)
    total = 0
    for sub, ext in [("ligands", "sdf"), ("pdbqt", "pdbqt"), ("poses", "pdbqt")]:
        fs = glob.glob(f"{D}/{sub}/*.{ext}")
        orph = [f for f in fs if os.path.basename(f)[: -(len(ext) + 1)] not in lib]
        total += len(orph)
        print(f"{sub:8s} {len(fs):4d} files, {len(orph):4d} not in library"
              f"{'' if dry_run else ' (deleted)'}")
        if not dry_run:
            for f in orph:
                os.remove(f)
    if dry_run and total:
        print(f"{total} orphan files; call prune_orphans(dry_run=False) to remove")
    return total
