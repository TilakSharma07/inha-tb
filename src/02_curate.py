"""Step 2 - filter cascade + molecule standardisation + replicate aggregation.

Every filter prints its survivor count, so the cascade in the README and in Figure 1a is
auditable end to end against a fresh run.

Two decisions worth stating explicitly:

  * Replicates are aggregated to the MEDIAN, not the mean. Whole-cell MICs are pooled from
    many laboratories and the distribution has outliers; the median is the robust summary.
  * The standard deviation across replicates is RETAINED (`pchembl_sd`). Its median over
    replicated molecules is the experimental noise floor - no model should be expected to
    predict below it, and it is the correct scale for reading model errors.
"""
import json, os
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
KEEP_TYPES = {"IC50", "Ki", "Kd", "MIC", "MIC90", "EC50"}
KEEP_UNITS = {"nM", "uM", "ug.mL-1"}
DATASETS = {"inha": "InhA_enzyme", "mtb": "Mtb_whole_cell"}   # keys in chembl_raw.json
KEEP_COLS = ["molecule_chembl_id", "canonical_smiles", "standard_type", "standard_value",
             "standard_units", "standard_relation", "pchembl_value", "assay_chembl_id",
             "assay_description", "data_validity_comment", "document_year",
             "document_chembl_id", "target_chembl_id", "potential_duplicate", "dataset",
             "molecule_pref_name"]


def curate():
    """Apply the filter cascade. Returns (dataframe, cascade) where cascade is the
    [(stage, surviving records)] list that Figure 1a and the README table are built from."""
    # the repo ships the gzipped pull (231 KB vs 6.6 MB) so the pipeline is reproducible
    # from step 2 without re-querying ChEMBL; a fresh run of step 1 writes the plain file
    gz = f"{D}/chembl_raw.json.gz"
    if os.path.exists(f"{D}/chembl_raw.json"):
        raw = json.load(open(f"{D}/chembl_raw.json"))
    elif os.path.exists(gz):
        import gzip
        with gzip.open(gz, "rt") as fh:
            raw = json.load(fh)
    else:
        raise SystemExit("no data/chembl_raw.json[.gz] - run src/01_fetch_chembl.py first")
    df = pd.DataFrame([dict(r, dataset=label)
                       for key, label in DATASETS.items() for r in raw[key]])
    df = df.reindex(columns=[c for c in KEEP_COLS if c in df.columns])
    cascade = []

    def stage(name, d):
        cascade.append((name, len(d)))
        print(f"{name:<34} {len(d):>6}")
        return d

    df = stage("ChEMBL records retrieved", df)
    df = stage("no data-validity flag", df[df.data_validity_comment.isna()])
    df = stage("exact relation (=)", df[df.standard_relation == "="])
    df = stage("potency endpoint", df[df.standard_type.isin(KEEP_TYPES)])
    df = stage("usable units", df[df.standard_units.isin(KEEP_UNITS)])
    df = stage("has pChEMBL value", df[df.pchembl_value.notna()])
    df["pchembl"] = df.pchembl_value.astype(float)
    df = stage("pChEMBL in [3,12]", df[df.pchembl.between(3, 12)])
    df = df[df.canonical_smiles.notna()]
    return df, cascade


def standardise(smi, lfc=rdMolStandardize.LargestFragmentChooser(),
                unc=rdMolStandardize.Uncharger()):
    """Largest fragment (drops salts/solvates), neutralise, canonical SMILES."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(unc.uncharge(lfc.choose(m)))


def main():
    df, cascade = curate()
    df["std_smiles"] = [standardise(s) for s in df.canonical_smiles]
    df = df[df.std_smiles.notna()]
    df.to_csv(f"{D}/chembl_filtered.csv", index=False)

    agg = (df.groupby(["dataset", "std_smiles"])
             .agg(molecule_chembl_id=("molecule_chembl_id", "first"),
                  pchembl=("pchembl", "median"),
                  pchembl_sd=("pchembl", "std"),
                  n_meas=("pchembl", "size"),
                  endpoints=("standard_type", lambda s: "|".join(sorted(set(s)))),
                  n_assay=("assay_chembl_id", "nunique"),
                  year_min=("document_year", "min"),
                  year_max=("document_year", "max"),
                  pref_name=("molecule_pref_name", "first"))
             .reset_index()
             [["dataset", "molecule_chembl_id", "std_smiles", "pchembl", "pchembl_sd",
               "n_meas", "endpoints", "n_assay", "year_min", "year_max", "pref_name"]])
    agg.to_csv(f"{D}/activity_aggregated.csv", index=False)

    cascade.append(("unique molecules\n(aggregated)", int(agg.std_smiles.nunique())))
    json.dump(dict(cascade=cascade, unit="unique standardised SMILES"),
              open(f"{D}/curation_cascade.json", "w"), indent=1)

    rep = agg[agg.n_meas > 1]
    print(f"\nunique molecules (std. SMILES)     {agg.std_smiles.nunique():>6}")
    print(f"molecule x dataset rows            {len(agg):>6}")
    print(f"replicated molecules               {len(rep):>6}")
    print(f"median inter-assay SD              {rep.pchembl_sd.median():>6.2f} log units "
          f"(experimental noise floor)")
    for ds, g in agg.groupby("dataset"):
        print(f"  {ds}: {len(g)} molecules, {int((g.pchembl >= 6).sum())} active at pChEMBL >= 6")


if __name__ == "__main__":
    main()
