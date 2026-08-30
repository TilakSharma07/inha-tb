"""Step 6 - select the docking library and generate 3D conformers.

Docking every curated molecule is unnecessary: what the screen has to answer is whether the
Vina score separates measured actives from measured inactives. That needs a subset that
(a) has measured potency, (b) spans both sides of the activity threshold, and (c) is
chemically diverse enough that the answer is not one congeneric series.

**Selection.** All InhA-enzyme molecules (the direct-target dataset, where a docking score
is mechanistically interpretable), plus a *class-balanced, scaffold-stratified* sample of
whole-cell molecules: N_PER_CLASS actives and N_PER_CLASS inactives, each drawn one-per-
scaffold first and then filled at random.

Balancing the two classes the *same* way matters. An earlier version sampled the whole-cell
arm without regard to class and the arm came out 100 % active, which makes discrimination
untestable on that dataset — ROC-AUC is undefined with one class. Sampling actives by
potency instead would be worse than useless: it would stack the most potent compounds on one
side and inflate the apparent separation.

**Ligand ids are content-addressed** (`L` + truncated hash of dataset|mol_id), not positional.
Positional ids meant any change to the selection renumbered every molecule and orphaned every
pose file already computed; hashed ids let a re-selection reuse the poses it overlaps with.

One conformer per molecule, ETKDG-embedded and MMFF-minimised. Conformer choice is a real
limitation: Vina treats the ligand as flexible about its rotatable bonds, but the starting
geometry still matters and a single conformer is the cheap option.
"""
import hashlib, os
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
N_PER_CLASS, SEED, CUTOFF = 150, 42, 6.0
# Vina's scoring function is parameterised on drug-like protein-ligand complexes, and at a
# fixed exhaustiveness the conformational search under-samples very flexible ligands, so a
# score for a 1,165 Da / 23-rotatable-bond molecule is not comparable to one for a 350 Da
# fragment. Molecules outside this envelope are excluded from the library rather than
# scored and quietly compared. The cut is on ligand properties only - never on the score -
# and is applied before docking, so it cannot select on the outcome.
MAX_MW, MAX_ROTB = 600.0, 12


def lig_id(dataset, mol_id):
    """Content-addressed ligand id: stable across re-selection, so poses survive."""
    h = hashlib.sha1(f"{dataset}|{mol_id}".encode()).hexdigest()[:8]
    return f"L{h}"


def embed(smi, seed=SEED):
    """SMILES -> 3D mol with hydrogens, ETKDG + MMFF. None on failure."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    if AllChem.EmbedMolecule(m, p) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    return m


def stratified(sub, n, seed=SEED):
    """One molecule per scaffold first (chemotype diversity), then fill at random."""
    first = sub.groupby("scaffold", sort=True).head(1)
    if len(first) >= n:
        return first.sample(n, random_state=seed)
    rest = sub.drop(first.index)
    return pd.concat([first, rest.sample(min(n - len(first), len(rest)), random_state=seed)])


def main():
    F = pd.read_csv(f"{D}/data/features.csv")
    n_all = len(F)
    F = F[(F.MW <= MAX_MW) & (F.RotB <= MAX_ROTB)].copy()
    print(f"dockable envelope (MW <= {MAX_MW:.0f}, RotB <= {MAX_ROTB}): "
          f"{len(F)}/{n_all} curated molecules")

    inha = F[F.dataset == "InhA_enzyme"].copy()
    wc = F[F.dataset == "Mtb_whole_cell"].copy()
    wc["active"] = (wc.pchembl >= CUTOFF).astype(int)

    pick = pd.concat([stratified(wc[wc.active == 1], N_PER_CLASS),
                      stratified(wc[wc.active == 0], N_PER_CLASS)])
    lib = pd.concat([inha, pick.drop(columns="active")]).reset_index(drop=True)
    lib["lig_id"] = [lig_id(d, m) for d, m in zip(lib.dataset, lib.mol_id)]
    assert lib.lig_id.is_unique, "ligand id collision"

    os.makedirs(f"{D}/dock/ligands", exist_ok=True)
    ok, fail, reused = [], [], 0
    for r in lib.itertuples():
        sdf = f"{D}/dock/ligands/{r.lig_id}.sdf"
        if os.path.exists(sdf) and os.path.getsize(sdf) > 0:
            ok.append(r.lig_id)
            reused += 1
            continue
        m = embed(r.smiles)
        if m is None:
            fail.append(r.lig_id)
            continue
        Chem.MolToMolFile(m, sdf)
        ok.append(r.lig_id)

    lib = lib[lib.lig_id.isin(ok)]
    lib.to_csv(f"{D}/dock/library.csv", index=False)
    act = (lib.pchembl >= CUTOFF)
    print(f"library: {len(lib)} molecules "
          f"({(lib.dataset=='InhA_enzyme').sum()} InhA enzyme, "
          f"{(lib.dataset=='Mtb_whole_cell').sum()} whole-cell) | "
          f"conformers reused {reused}, embed failures {len(fail)}")
    print(f"actives at pChEMBL >= {CUTOFF}: {int(act.sum())} / {len(lib)} "
          f"({100*act.mean():.0f} %)")
    for ds, sub in lib.groupby("dataset"):
        a = (sub.pchembl >= CUTOFF).sum()
        print(f"  {ds}: {len(sub)} molecules, {a} active ({100*a/len(sub):.0f} %), "
              f"{sub.scaffold.nunique()} scaffolds")


if __name__ == "__main__":
    main()
