"""Step 3 - Morgan fingerprints + physicochemical descriptors + Bemis-Murcko scaffold split.

Split convention (the part that is easy to get backwards): the LARGEST scaffold groups go
to TRAIN and the rare/singleton scaffolds go to TEST. This is the DeepChem convention and
the harder, realistic setting - the model must generalise to chemotypes it has never seen.

Inverting it puts the big congeneric series in the test set, which then scores HIGHER than
a random split. If you ever see scaffold > random, the split is backwards, not the model
generalising.
"""
import json, os
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
NBITS, RADIUS, TEST_FRAC, CUTOFF, SEED = 2048, 2, 0.25, 6.0, 42

# The 17 descriptors concatenated onto the fingerprint block.
DCOLS = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "RingCount", "AromRings",
         "FracCsp3", "HeavyAtoms", "FormalCharge", "QED", "BertzCT",
         "NumHeteroatoms", "StereoCentres", "Lipinski_viol", "Veber_ok"]


def descriptors(m):
    return dict(
        MW=Descriptors.MolWt(m), cLogP=Crippen.MolLogP(m),
        TPSA=rdMolDescriptors.CalcTPSA(m), HBD=rdMolDescriptors.CalcNumHBD(m),
        HBA=rdMolDescriptors.CalcNumHBA(m),
        RotB=rdMolDescriptors.CalcNumRotatableBonds(m),
        RingCount=rdMolDescriptors.CalcNumRings(m),
        AromRings=rdMolDescriptors.CalcNumAromaticRings(m),
        FracCsp3=rdMolDescriptors.CalcFractionCSP3(m),
        HeavyAtoms=m.GetNumHeavyAtoms(),
        FormalCharge=Chem.GetFormalCharge(m), QED=QED.qed(m),
        BertzCT=Descriptors.BertzCT(m),
        NumHeteroatoms=rdMolDescriptors.CalcNumHeteroatoms(m),
        StereoCentres=rdMolDescriptors.CalcNumAtomStereoCenters(m),
    )


def scaffold_split(scaffolds, test_frac=TEST_FRAC):
    """Largest scaffold groups -> train, rare scaffolds -> test. Returns (train, test) index lists.

    Deterministic: groups are ordered by descending size with the scaffold SMILES as
    tie-breaker, so the result does not depend on dict iteration order.
    """
    groups = {}
    for i, s in enumerate(scaffolds):
        groups.setdefault(s, []).append(i)
    order = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    n_test = int(len(scaffolds) * test_frac)
    train, test = [], []
    for _, idx in order:                       # biggest first -> train until train quota met
        if len(train) < len(scaffolds) - n_test:
            train += idx
        else:
            test += idx
    return sorted(train), sorted(test)


def main():
    agg = pd.read_csv(f"{D}/activity_aggregated.csv")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)

    recs, fps = [], []
    for r in agg.itertuples():
        m = Chem.MolFromSmiles(r.std_smiles)
        if m is None:
            continue
        d = descriptors(m)
        d.update(dataset=r.dataset, mol_id=r.molecule_chembl_id, smiles=r.std_smiles,
                 pchembl=r.pchembl, n_meas=r.n_meas, endpoints=r.endpoints,
                 scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=m))
        recs.append(d)
        fps.append(np.array(gen.GetFingerprint(m), dtype=np.uint8))

    F = pd.DataFrame(recs)
    F["Lipinski_viol"] = ((F.MW > 500).astype(int) + (F.cLogP > 5).astype(int)
                          + (F.HBD > 5).astype(int) + (F.HBA > 10).astype(int))
    F["Veber_ok"] = ((F.RotB <= 10) & (F.TPSA <= 140)).astype(int)
    X = np.vstack(fps)
    assert X.shape == (len(F), NBITS), (X.shape, len(F))

    F.to_csv(f"{D}/features.csv", index=False)
    np.save(f"{D}/fp_matrix.npy", X)
    print(f"featurised {len(F)} molecule x dataset rows | fp {X.shape} | "
          f"mean bits set {X.sum(1).mean():.0f}")

    splits = {}
    for ds, sub in F.groupby("dataset"):
        gi = sub.index.to_numpy()                       # positions into F / fp_matrix
        y = (sub.pchembl >= CUTOFF).astype(int).to_numpy()
        scaf = np.array(sub.scaffold.tolist())
        tr, te = scaffold_split(scaf.tolist())
        rng = np.random.default_rng(SEED)
        perm = rng.permutation(len(gi))
        rte, rtr = sorted(perm[:len(te)].tolist()), sorted(perm[len(te):].tolist())
        leak_s = len(set(scaf[tr]) & set(scaf[te]))
        leak_r = len(set(scaf[rtr]) & set(scaf[rte]))
        splits[ds] = dict(global_idx=gi.tolist(), y=y.tolist(),
                          scaf_train=tr, scaf_test=te, rand_train=rtr, rand_test=rte,
                          leaked_scaffolds_scaffold=leak_s, leaked_scaffolds_random=leak_r,
                          n_scaffolds=int(len(set(scaf))))
        print(f"{ds}: n={len(gi)} actives={int(y.sum())} scaffolds={len(set(scaf))} | "
              f"scaffold test={len(te)} (leaked {leak_s}) | random test={len(rte)} (leaked {leak_r})")

    json.dump(dict(cutoff=CUTOFF, nbits=NBITS, radius=RADIUS, test_frac=TEST_FRAC,
                   descriptor_cols=DCOLS, splits=splits),
              open(f"{D}/splits.json", "w"))


if __name__ == "__main__":
    main()
