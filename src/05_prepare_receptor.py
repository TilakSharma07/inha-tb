"""Step 5 - build the docking receptor and box from PDB 4TZK.

Two things here are easy to get wrong and are handled explicitly:

1. NCBI MMDB-served PDB files interleave biological-assembly copies as MODEL blocks, and
   the HETATM records sit AFTER `ENDMDL` on non-standard chain ids ('N' for NAD, '(' for
   the inhibitor). A naive MODEL/ENDMDL parse silently drops the cofactor and ligand.
2. The NAD cofactor is RETAINED in the receptor. InhA is NADH-dependent and the inhibitor
   stacks on the nicotinamide ring - stripping HETATMs wholesale destroys the pocket.
"""
import json, os
import numpy as np

D = os.path.join(os.path.dirname(__file__), "..")
PDB_ID, PAD, MIN_BOX = "4TZK", 5.0, 22.0


def xyz(rows):
    return np.array([[float(r[30:38]), float(r[38:46]), float(r[46:54])] for r in rows])


def main():
    L = open(f"{D}/structures/{PDB_ID}.pdb").read().split("\n")
    mdl = [i for i, l in enumerate(L) if l.startswith("MODEL")]
    end = [i for i, l in enumerate(L) if l.startswith("ENDMDL")]
    # assembly copy 1 = protein between MODEL[0]..ENDMDL[0], HETATM block up to MODEL[1]
    prot = [l for l in L[mdl[0] + 1:end[0]] if l.startswith("ATOM")]
    blk = [l for l in L[end[0]:mdl[1]] if l.startswith("HETATM")]
    nad = [l for l in blk if l[17:20].strip() == "NAD"]
    lig = [l for l in blk if l[17:20].strip() == "641"]
    print(f"copy 1: protein={len(prot)}  NAD={len(nad)}  ligand={len(lig)}")
    assert prot and nad and lig, "cofactor or ligand not found - check the MODEL layout"

    os.makedirs(f"{D}/dock", exist_ok=True)
    nad_A = [l[:21] + "A" + l[22:] for l in nad]      # unify chain id
    open(f"{D}/dock/receptor_raw.pdb", "w").write("\n".join(prot + nad_A) + "\nEND\n")
    open(f"{D}/dock/ref_ligand.pdb", "w").write("\n".join(lig) + "\nEND\n")

    lx, nx, rx = xyz(lig), xyz(nad_A), xyz(prot + nad_A)
    ctr, span = lx.mean(0), lx.max(0) - lx.min(0)
    box = np.maximum(span + 2 * PAD, MIN_BOX).round(1)
    inside = int(np.all(np.abs(rx - ctr) <= box / 2, axis=1).sum())
    dmin = float(np.linalg.norm(lx[:, None] - nx[None], axis=2).min())
    print(f"box {box.tolist()} A   receptor atoms enclosed={inside}   min ligand-NAD={dmin:.2f} A")

    json.dump(dict(center=ctr.round(3).tolist(), size=box.tolist(), receptor_pdb=PDB_ID,
                   resolution_A=1.62, ref_ligand_resname="641", atoms_in_box=inside,
                   min_lig_nad_A=round(dmin, 2),
                   cofactor="NAD retained - InhA is NADH-dependent"),
              open(f"{D}/dock/box.json", "w"), indent=1)

    os.system(f'obabel -ipdb {D}/dock/receptor_raw.pdb -opdbqt -O {D}/dock/receptor.pdbqt '
              f'-xr -p 7.4 --partialcharge gasteiger')


if __name__ == "__main__":
    main()
