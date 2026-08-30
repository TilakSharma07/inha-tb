"""Dock a prepared ligand library into InhA (PDB 4TZK) with AutoDock Vina.
Parallel: N worker processes, each running Vina with a small cpu budget."""
import glob, json, os, sys, time
import multiprocessing as mp
import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
BOX = json.load(open(f"{D}/box.json"))
NPROC, VINA_CPU, EXH = 6, 2, 12

def prep_one(sdf):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    lid = os.path.basename(sdf)[:-4]
    out = f"{D}/pdbqt/{lid}.pdbqt"
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return lid
    m = Chem.MolFromMolFile(sdf, removeHs=False)
    if m is None: return None
    try:
        setups = MoleculePreparation().prepare(m)
        txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok: return None
        open(out, "w").write(txt)
        return lid
    except Exception:
        return None

def dock_chunk(args):
    idx, chunk = args
    from vina import Vina
    v = Vina(sf_name="vina", cpu=VINA_CPU, seed=42, verbosity=0)
    v.set_receptor(f"{D}/receptor.pdbqt")
    v.compute_vina_maps(center=BOX["center"], box_size=BOX["size"])
    out = []
    for j, lid in enumerate(chunk):
        try:
            v.set_ligand_from_file(f"{D}/pdbqt/{lid}.pdbqt")
            v.dock(exhaustiveness=EXH, n_poses=9)
            en = v.energies(n_poses=9)
            v.write_poses(f"{D}/poses/{lid}.pdbqt", n_poses=3, overwrite=True)
            out.append(dict(lig_id=lid, vina_kcal=float(en[0][0]),
                            vina_mean_top3=float(np.mean([e[0] for e in en[:3]])),
                            pose_spread=float(np.std([e[0] for e in en])), n_poses=len(en)))
        except Exception as e:
            out.append(dict(lig_id=lid, vina_kcal=None, error=type(e).__name__))
        if idx == 0 and (j+1) % 10 == 0:
            print(f"  worker0 {j+1}/{len(chunk)}", flush=True)
    return out

if __name__ == "__main__":
    os.makedirs(f"{D}/pdbqt", exist_ok=True); os.makedirs(f"{D}/poses", exist_ok=True)
    sdfs = sorted(glob.glob(f"{D}/ligands/*.sdf"))
    t0 = time.time()
    with mp.Pool(NPROC) as pool:
        ids = [x for x in pool.map(prep_one, sdfs) if x]
    print(f"prepared {len(ids)}/{len(sdfs)} pdbqt in {time.time()-t0:.0f}s", flush=True)

    chunks = [(i, list(c)) for i, c in enumerate(np.array_split(np.array(ids), NPROC))]
    t1 = time.time()
    with mp.Pool(NPROC) as pool:
        rows = [r for part in pool.map(dock_chunk, chunks) for r in part]
    json.dump(rows, open(f"{D}/docking_scores.json", "w"))
    good = [r["vina_kcal"] for r in rows if r.get("vina_kcal") is not None]
    print(f"DONE {len(good)}/{len(ids)} docked in {(time.time()-t1)/60:.1f} min", flush=True)
    print(f"vina: min={min(good):.2f} median={np.median(good):.2f} max={max(good):.2f} kcal/mol", flush=True)
