"""Step 6 - the control that decides whether the docking setup can be trusted.

Redock the co-crystallised ligand and measure RMSD to its experimental pose. Anything
above ~2 A means the box, protonation, or receptor prep is wrong and every downstream
score is noise.

RMSD note: the docked pose and the crystal pose do NOT share atom ordering, and partial
sanitisation of a PDBQT-derived molecule leaves implicit-H counts at zero, so RDKit's
template/substructure routes fail. The measurement used here is an element-constrained
optimal assignment (Hungarian) over heavy atoms in absolute coordinates - no realignment.
"""
import glob, json, os, subprocess
import numpy as np
from scipy.optimize import linear_sum_assignment
from rdkit import Chem, RDLogger
from vina import Vina
RDLogger.DisableLog("rdApp.*")

D = os.path.join(os.path.dirname(__file__), "..", "dock")


def heavy_from_sdf(f):
    L = open(f).read().split("\n")
    na = int(L[3][:3])
    rows = [(l[31:34].strip(), float(l[0:10]), float(l[10:20]), float(l[20:30])) for l in L[4:4 + na]]
    rows = [r for r in rows if r[0] != "H"]
    return np.array([r[0] for r in rows]), np.array([[r[1], r[2], r[3]] for r in rows])


def rmsd_assign(el_p, P, el_r, Rx):
    Dm = np.linalg.norm(P[:, None, :] - Rx[None, :, :], axis=2)
    Dm = np.where(el_p[:, None] == el_r[None, :], Dm, 1e6)
    i, j = linear_sum_assignment(Dm)
    d = Dm[i, j]
    return float(np.sqrt((d ** 2).mean())) if d.max() < 1e5 else np.nan


def main():
    box = json.load(open(f"{D}/box.json"))
    subprocess.run(["obabel", "-ipdb", f"{D}/ref_ligand.pdb", "-osdf",
                    "-O", f"{D}/ref_ligand.sdf", "-p", "7.4"], capture_output=True)
    subprocess.run(["mk_prepare_ligand.py", "-i", f"{D}/ref_ligand.sdf",
                    "-o", f"{D}/ref_ligand.pdbqt"], capture_output=True)

    v = Vina(sf_name="vina", cpu=8, seed=42, verbosity=0)
    v.set_receptor(f"{D}/receptor.pdbqt")
    v.set_ligand_from_file(f"{D}/ref_ligand.pdbqt")
    v.compute_vina_maps(center=box["center"], box_size=box["size"])
    crystal_score = float(v.score()[0])
    v.dock(exhaustiveness=32, n_poses=20)
    v.write_poses(f"{D}/ref_redock.pdbqt", n_poses=20, overwrite=True)
    en = v.energies()

    for f in glob.glob(f"{D}/ref_redock[0-9]*.sdf"):
        os.remove(f)
    subprocess.run(["obabel", "-ipdbqt", f"{D}/ref_redock.pdbqt", "-osdf",
                    "-O", f"{D}/ref_redock.sdf", "-m"], capture_output=True)

    ref = Chem.MolFromMolFile(f"{D}/ref_ligand.sdf", removeHs=True)
    el_r = np.array([a.GetSymbol() for a in ref.GetAtoms()])
    Rx = ref.GetConformer().GetPositions()

    files = sorted(glob.glob(f"{D}/ref_redock[0-9]*.sdf"),
                   key=lambda f: int("".join(c for c in os.path.basename(f) if c.isdigit())))
    rms = []
    for f in files:
        el_p, P = heavy_from_sdf(f)
        rms.append(rmsd_assign(el_p, P, el_r, Rx) if len(P) == len(Rx) else np.nan)

    ok = [r for r in rms if r == r]
    print(f"crystal-pose score {crystal_score:.2f} kcal/mol | redock top {en[0][0]:.2f}")
    print(f"RMSD top-ranked pose {rms[0]:.2f} A | best {min(ok):.2f} A | <2A: {sum(r<2 for r in ok)}/{len(ok)}")
    json.dump(dict(crystal_score=crystal_score, redock_top5=[float(e[0]) for e in en[:5]],
                   rmsd_top=rms[0], rmsd_best=min(ok), rmsd_all=[round(float(r), 3) for r in rms],
                   method="element-constrained optimal assignment, heavy atoms, absolute coords"),
              open(f"{D}/redock_validation.json", "w"), indent=1)
    assert rms[0] < 2.0, "REDOCK FAILED - do not trust downstream scores"


if __name__ == "__main__":
    main()
