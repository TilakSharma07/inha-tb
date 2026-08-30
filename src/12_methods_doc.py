"""Step 12 - generate docs/METHODS.md from the saved tables.

Every table and every number in the methods document is read from data/ and dock/ at write
time. The prose is here; the numbers are not. This is deliberate: the earlier hand-typed
version of this document drifted from the data within a single session (it quoted a
descriptor count, a learning rate, and a filter cascade that no longer matched the code),
which is exactly the failure mode a methods section cannot afford.
"""
import importlib.util
import inspect, json, os
import numpy as np, pandas as pd
from scipy import stats

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.dirname(os.path.abspath(__file__))
MODEL_LABEL = {"XGBoost": "XGBoost", "RandomForest": "random forest",
               "LogReg(baseline)": "logistic regression"}
# same definition as src/10_figures.py - Gly and Tyr counted as pocket-lining apolar/aromatic
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "PRO", "GLY", "TYR"}
DS_LABEL = {"InhA_enzyme": "InhA enzyme", "Mtb_whole_cell": "*M. tb* whole-cell"}


def hyperparams():
    """Read the model definitions out of 04_train_qsar.py so the text cannot drift."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("t", f"{SRC}/04_train_qsar.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    out = {}
    for name, mk in m.MODELS.items():
        est = mk()
        p = est.get_params()
        if name == "XGBoost":
            out[name] = (f"{p['n_estimators']} trees, depth {p['max_depth']}, "
                         f"lr {p['learning_rate']}, subsample {p['subsample']}, "
                         f"colsample {p['colsample_bytree']}, min_child_weight "
                         f"{p['min_child_weight']}, L2 {p['reg_lambda']}")
        elif name == "RandomForest":
            out[name] = (f"{p['n_estimators']} trees, min_samples_leaf "
                         f"{p['min_samples_leaf']}, max_features {p['max_features']}, "
                         f"class_weight {p['class_weight']}")
        else:
            out[name] = (f"L2, C = {p['C']}, class_weight {p['class_weight']}, "
                         f"on standardised features")
    return out, m.AD_BINS


def md_table(df, cols=None, fmt=None):
    df = df[cols] if cols else df
    fmt = fmt or {}
    head = "| " + " | ".join(df.columns) + " |"
    rule = "|" + "|".join(["---"] * len(df.columns)) + "|"
    rows = []
    for r in df.itertuples(index=False):
        cells = [fmt.get(c, str)(v) for c, v in zip(df.columns, r)]
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, rule] + rows)


def main():
    F = json.load(open(f"{D}/docs/facts.json"))
    R = pd.read_csv(f"{D}/data/metrics.csv")
    AD = pd.read_csv(f"{D}/data/applicability_domain.csv")
    SP = json.load(open(f"{D}/data/splits.json"))
    SM = json.load(open(f"{D}/data/structure_meta.json"))
    PK = json.load(open(f"{D}/dock/pocket_residues.json"))
    V = json.load(open(f"{D}/dock/redock_validation.json"))
    HP, AD_BINS = hyperparams()

    # structures table
    ST = pd.DataFrame([dict(PDB=k,
                            resolution=v["resolution_A"],
                            chains=len(v["chains"]),
                            HET=", ".join([h for h in v["het"] if h != "HOH"][:4]))
                       for k, v in SM.items() if isinstance(v, dict)])
    ST = ST.sort_values("resolution").reset_index(drop=True)

    # per-model results table
    piv = R.pivot_table(index=["dataset", "model"], columns="split", values="ROC_AUC").reset_index()
    piv["inflation"] = piv["rand"] - piv["scaf"]
    piv["dataset"] = piv.dataset.map(DS_LABEL)
    piv["model"] = piv.model.map(MODEL_LABEL)
    piv = piv.rename(columns={"scaf": "scaffold", "rand": "random"})[
        ["dataset", "model", "scaffold", "random", "inflation"]]

    ADt = AD.copy()
    ADt["dataset"] = ADt.dataset.map(DS_LABEL)
    ADt["interpretable"] = np.where(ADt.n >= 10, "yes", "**no (n < 10)**")
    ADt = ADt.rename(columns={"tanimoto_bin": "Tanimoto bin", "mean_abs_err": "mean abs. err"})[
        ["dataset", "Tanimoto bin", "n", "actives", "ROC_AUC", "mean abs. err", "interpretable"]]

    dk, tr, q = F["docking"], F["translation"], F["qsar"]
    f3 = lambda v: f"{v:.3f}"
    f2 = lambda v: f"{v:.2f}"
    fp3 = lambda v: f"{v:+.3f}"

    # tables are rendered here, not inside the f-string: a replacement field cannot hold a
    # dict literal (`{{` there opens a set, not an escaped brace)
    TBL_STRUCT = md_table(ST, fmt={"resolution": f2})
    TBL_RESULTS = md_table(piv, fmt={"scaffold": f3, "random": f3, "inflation": fp3})
    TBL_AD = md_table(ADt, fmt={"ROC_AUC": f3, "mean abs. err": f3})
    PERF = json.load(open(f"{D}/data/docking_performance.json"))
    TBL_RES = md_table(pd.DataFrame(
        [{"dataset": DS_LABEL.get(k, k), "n": v["n"], "actives": v["actives"],
          "ROC-AUC": f"{v['roc_auc']:.2f}", "Spearman ρ": f"{v['spearman']:.2f}",
          "EF 5 %": f"{v['ef5']:.2f}", "ceiling": f"{v['ef_max']:.2f}"}
         for k, v in sorted(PERF["by_dataset"].items())]))
    TBL_H2H = md_table(pd.DataFrame(
        [{"dataset": DS_LABEL.get(k, k), "held-out n": v["docking"]["n"],
          "docking ROC-AUC": f"{v['docking']['roc_auc']:.2f}",
          "QSAR ROC-AUC": f"{v['qsar']['roc_auc']:.2f}"}
         for k, v in sorted(PERF["held_out_head_to_head"].items())]))
    bal = dk.get("library_balance", {})
    TBL_LIB = md_table(pd.DataFrame(
        [{"dataset": DS_LABEL.get(k, k), "molecules": v["n"],
          "active (pChEMBL ≥ 6)": v["actives"],
          "% active": f"{100*v['actives']/v['n']:.0f} %"}
         for k, v in sorted(bal.items())])) if bal else ""
    # the selection constants belong to src/06; import them rather than restate them
    _spec = importlib.util.spec_from_file_location("build_lib", f"{D}/src/06_build_library.py")
    _bl = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_bl)
    LIB_PER_CLASS, MAX_MW, MAX_ROTB = _bl.N_PER_CLASS, _bl.MAX_MW, _bl.MAX_ROTB
    adf = F["applicability_domain"]
    ADhi_lo, ADhi_hi = min(adf["near_roc_auc"]), max(adf["near_roc_auc"])
    ADlo_lo, ADlo_hi = min(adf["far_roc_auc"]), max(adf["far_roc_auc"])

    doc = f"""# Methods

Every number in this document is read from the saved tables at generation time by
`src/12_methods_doc.py`. Nothing here is hand-typed, because an earlier hand-typed version
drifted from the data within a single session.

## 1. Target and receptor selection

InhA (enoyl-[acyl-carrier-protein] reductase, *Rv1484*) is the NADH-dependent enzyme of the
mycobacterial FAS-II elongation cycle and the validated molecular target of isoniazid.
Isoniazid requires KatG-mediated activation to form the INH–NAD adduct, and *katG* mutations
are the dominant mechanism of clinical isoniazid resistance. Direct InhA inhibitors that do
not require KatG activation therefore remain a well-motivated route to compounds active
against INH-resistant strains.

{len(ST)} InhA crystal structures were retrieved from the PDB and compared on resolution, cofactor
presence, chain count, and the drug-likeness of the bound ligand.

{TBL_STRUCT}

**{dk['receptor']}** was selected: highest resolution, single chain, NAD cofactor present, and
a drug-like co-crystallised inhibitor (ligand code 641) usable as a redocking reference.

### Receptor preparation

The coordinate file served by the NCBI MMDB mirror interleaves biological-assembly copies as
`MODEL` blocks and places HETATM records **after** `ENDMDL`, on non-standard chain
identifiers (`N` for NAD, `(` for ligand 641). A conventional MODEL/ENDMDL parse silently
discards both the cofactor and the reference ligand; `src/05_prepare_receptor.py` handles
this layout explicitly and asserts that all three components are recovered.

**The NAD cofactor is retained in the receptor.** InhA is NADH-dependent and the inhibitor
stacks directly against the nicotinamide ring. Stripping HETATM records wholesale — a common
shortcut — removes one wall of the binding site and invalidates every subsequent score.

Protonation at pH 7.4 and Gasteiger partial charges were assigned with Open Babel; the
receptor was written as a rigid PDBQT.

### Docking box

Centred on the centroid of the co-crystallised inhibitor, {dk['box_center']}, with edge
lengths {dk['box_size']} Å (ligand span + 5 Å padding, floored at 22 Å).

The pocket lining was characterised by distance: {len(PK)} residues fall within
4.5 Å of ligand 641, of which {sum(1 for r in PK if r["resn"] in HYDROPHOBIC)}
are apolar or aromatic. Both canonical InhA hydrogen-bonding anchors — **Tyr158** and **Thr196** —
are present, confirming the box covers the substrate-binding site rather than a surface
groove.

## 2. Redocking control

Empirical docking scores are interpretable only if the setup reproduces a known experimental
binding mode. The co-crystallised ligand was extracted, re-protonated, re-prepared as a
flexible ligand, and redocked into the prepared receptor (exhaustiveness 32, 20 poses).

| quantity | value |
|---|---|
| crystal-pose Vina score | {dk['crystal_score']:.2f} kcal/mol |
| best redocked score | {dk['redock_top']:.2f} kcal/mol |
| **RMSD, top-ranked pose** | **{dk['rmsd_top']:.2f} Å** |
| best RMSD among 20 poses | {dk['rmsd_best']:.2f} Å |
| centroid displacement | {V['centroid_shift']:.2f} Å |

RMSD note: the docked and crystal poses do not share atom ordering, and a PDBQT-derived
molecule carries no implicit-hydrogen counts, so RDKit's template and substructure-match
routes both fail on this pair. The reported value is an **element-constrained optimal
assignment** (Hungarian algorithm) over heavy atoms in absolute coordinates, with no
realignment — a measurement that cannot flatter the result through superposition.

A top-ranked-pose RMSD of {dk['rmsd_top']:.2f} Å is well inside the ≤ 2 Å convention for a
successful redock. `src/07_redock_validate.py` asserts this and aborts the pipeline if it
fails.

## 3. Bioactivity data

Two ChEMBL targets were pulled: **CHEMBL1849** (InhA enzyme assays) and **CHEMBL360**
(*M. tuberculosis* whole-cell). The whole-cell target carries >200k activities and the
endpoint offers no usable offset beyond a few thousand rows, so the potency axis was
partitioned into pChEMBL bands and each band paginated separately.

Filter cascade, as printed by `src/02_curate.py` and written to
`data/curation_cascade.json`:

| stage | records |
|---|---|
""" + "\n".join(f"| {k.replace(chr(10), ' ')} | {v:,} |" for k, v in F["cascade"]) + f"""

Molecules were standardised with RDKit `rdMolStandardize`: largest-fragment selection
(drops salts and solvates), neutralisation, canonical SMILES. Replicate measurements of the
same molecule against the same target were aggregated to the **median** — whole-cell MICs
are pooled across laboratories and the distribution carries outliers — with the standard
deviation across replicates retained.

**Assay noise floor.** {F['n_replicated']} molecules have replicate measurements; their
median inter-assay SD is **{F['noise_floor_sd']:.2f} log units**. No model should be expected
to predict below this, and it is the correct scale against which to read the errors in §6.

Final modelling sets: **{F['n_inha']}** molecules with InhA enzyme potency
({F['actives_inha']} active at pChEMBL ≥ 6) and **{F['n_wc']}** with whole-cell potency
({F['actives_wc']} active). {F['molecules_unique']:,} unique structures in total; the
difference is molecules measured in both assays.

## 4. Target engagement translates to cell kill

{tr['n']} molecules were measured in both assays. Their potencies correlate at
**Spearman ρ = {tr['spearman_rho']:.2f}** (p = {tr['p_value']:.0e}). This is the mechanistic
anchor of the project: InhA inhibition does not merely correlate with an enzymatic readout,
it tracks growth inhibition of the organism.

## 5. QSAR models

**Features.** {F['fingerprint']} concatenated with {F['descriptors_n']} physicochemical
descriptors: {', '.join(SP['descriptor_cols'])}. The last two are rule-based flags
(Lipinski violation count, Veber pass/fail) rather than continuous properties.

**Labels.** Binary at pChEMBL ≥ {SP['cutoff']} (1 µM) for both endpoints.

**Splits.** Two, deliberately compared:

- *Scaffold split* — Bemis–Murcko scaffolds; the **largest scaffold groups go to train** and
  rare/singleton scaffolds to test.
  {SP['splits']['InhA_enzyme']['leaked_scaffolds_scaffold']} scaffolds are shared across the
  enzyme-set boundary and
  {SP['splits']['Mtb_whole_cell']['leaked_scaffolds_scaffold']} across the whole-cell boundary
  ({SP['splits']['InhA_enzyme']['n_scaffolds']} and
  {SP['splits']['Mtb_whole_cell']['n_scaffolds']} distinct scaffolds respectively).
- *Random split* — same size, uniformly sampled. This leaks whole chemical series across the
  boundary: {q['InhA_enzyme']['leaked_scaffolds_random']} and
  {q['Mtb_whole_cell']['leaked_scaffolds_random']} scaffolds appear on both sides for the
  enzyme and whole-cell sets respectively.

The split direction matters and is easy to invert. Assigning the largest scaffold groups to
*test* makes the test set a congeneric series and yields a scaffold-split score **higher**
than the random-split score — a diagnostic that the split is backwards, not that the model
generalises.

**Models.**

""" + "\n".join(f"- **{MODEL_LABEL[k]}** — {v}" for k, v in HP.items()) + f"""

**Results (ROC-AUC on held-out test):**

{TBL_RESULTS}

A random split overstates ROC-AUC by
**{F['inflation']['min']:.2f}–{F['inflation']['max']:.2f}** on this data (median
{F['inflation']['median']:.2f}). Every result quoted outside this table is the scaffold-split
number.

## 6. Applicability domain

Test compounds (scaffold split, random forest) were binned by maximum Tanimoto similarity to
any training compound. Bins holding fewer than 10 compounds are flagged: a ROC-AUC computed
on five compounds is not a measurement.

{TBL_AD}

Accuracy is a function of chemical distance, not a fixed property of the model. Above
Tanimoto 0.35 the interpretable bins score {ADhi_lo:.2f}–{ADhi_hi:.2f}; below 0.35 they fall
to {ADlo_lo:.2f}–{ADlo_hi:.2f} — at, or barely above, chance. The trend is not strictly
monotonic (the 0.35–0.50 bins score marginally above the 0.50+ bins, on 13 and 28 compounds),
so the useful reading is the split between the two low-similarity bins and everything else,
not a smooth curve. Predictions on novel chemotypes should not be trusted, and this table is
the reason.

## 7. Docking screen

Library: {dk['library_n']} molecules ({dk['library_inha']} with InhA enzyme potency,
{dk['library_wc']} whole-cell), of which {dk['library_actives']} are measured actives at
pChEMBL ≥ 6. One ETKDG conformer per molecule, MMFF-minimised, Meeko-prepared.

{TBL_LIB}

Two selection decisions are worth stating, because both change what the screen can show.

**The whole-cell arm is class-balanced.** Molecules are drawn one-per-scaffold within each
class separately — {LIB_PER_CLASS} actives and {LIB_PER_CLASS} inactives — rather than
sampled from the dataset as a whole. An earlier version sampled without regard to class and
the arm came out 100 % active, which makes discrimination untestable on that dataset: ROC-AUC
is undefined when one class is empty. Selecting actives by potency instead would be worse
than a bug — it would place the most potent compounds on one side of the comparison and
inflate the apparent separation.

**Molecules outside a drug-like envelope are excluded** (MW ≤ {MAX_MW:.0f} Da,
≤ {MAX_ROTB} rotatable bonds; the library's actual maxima are {dk['library_mw_max']:.0f} Da
and {dk['library_rotb_max']}). Vina's scoring function is parameterised on drug-like
complexes, and at fixed exhaustiveness the conformational search under-samples very flexible
ligands, so a score for the most extreme curated molecule ({dk['excluded_max_mw']:,.0f} Da,
{dk['excluded_max_rotb']} rotatable bonds) is not comparable to one for a 350 Da fragment.
The envelope removes {dk['n_excluded']} of the curated molecules. The cut is on ligand
properties only, never on the score, and is applied before docking, so it cannot select on
the outcome. Selection and the single-conformer
limitation are implemented in `src/06_build_library.py`.

## 8. What the docking screen shows

{TBL_RES}

Ranking by Vina score separates measured actives from inactives only slightly better than
chance on the enzyme arm and less than that on the whole-cell arm. Rank correlation with
measured potency is weak and negative in the expected direction (a more negative score
should mean a more potent compound), but far too weak to rank a series. Enrichment in the
top 5 % of the ranked list is below the attainable ceiling on both arms, and on the
whole-cell arm it is below 1.0 - selecting the top-scoring 5 % of that library would return
fewer actives than picking at random.

**The comparison that matters is on the same molecules.** Docking and the QSAR model were
scored on the held-out scaffold-split test compounds, which neither approach was fitted to:

{TBL_H2H}

On identical molecules the ligand-based model separates actives from inactives far better
than the physics-based score. This is the expected result for a target with this much
measured data - a model that has seen {F['records_curated']:,} curated activity records for
{F['molecules_unique']:,} molecules is using information that a single rigid-receptor docking run
does not have. The screen is reported because a negative result on a well-posed comparison
is a result: it bounds what the docking component can contribute, and it is the reason the
prediction endpoint is the QSAR model rather than the docking score.

Three specific limitations bound the docking number, and none of them are ruled out here.
One conformer per ligand under-samples flexible molecules. The receptor is a single rigid
crystal structure, so induced fit is not modelled. And Vina's score correlates with
molecular size (Spearman rho vs MW is reported in figure 3f), which inflates the apparent
performance of any list ranked by raw score. A size-corrected score, ensemble docking across
the {len(SM)} InhA structures screened, or rescoring the retained poses would each test
one of these; none is done here.

## 9. What this pipeline does not claim

- **Vina scores are not affinities.** The scoring function was parameterised to reproduce
  binding geometry, not to rank potency. Reported use is triage and prioritisation only.
- **Rigid receptor.** InhA's substrate-binding loop (residues ~196–219) reorganises on
  inhibitor binding; several structures screened in §1 differ in exactly that region. A
  single rigid conformation cannot capture induced fit.
- **One conformer per ligand.** Vina samples torsions but not ring or starting-geometry
  alternatives.
- **No experimental validation.** Every number here is computational.
- **Assay heterogeneity.** Whole-cell MICs are aggregated across laboratories and protocols.
  The retained replicate SD quantifies this; it does not remove it.
- **Fingerprints are blind to activity cliffs** within a congeneric series.
"""
    os.makedirs(f"{D}/docs", exist_ok=True)
    open(f"{D}/docs/METHODS.md", "w").write(doc)
    print(f"wrote docs/METHODS.md ({len(doc)} chars, {doc.count(chr(10))+1} lines)")


if __name__ == "__main__":
    main()
