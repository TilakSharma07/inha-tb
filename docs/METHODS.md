# Methods

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

10 InhA crystal structures were retrieved from the PDB and compared on resolution, cofactor
presence, chain count, and the drug-likeness of the bound ligand.

| PDB | resolution | chains | HET |
|---|---|---|---|
| 4TZK | 1.62 | 1 | NAD, 641 |
| 4OIM | 1.85 | 1 | NAD, ACT, JUS |
| 2NSD | 1.90 | 1 | NAD, 4PI |
| 2AQ8 | 1.92 | 1 | NAI, LYS |
| 3FNG | 1.97 | 1 | NAD, JPL |
| 5COQ | 2.30 | 4 | NAD, TCU, NA |
| 2B35 | 2.30 | 2 | NAD, TCL |
| 4OXY | 2.35 | 4 | NAD, 1TN |
| 4COD | 2.40 | 4 | NAD, KV1 |
| 1ZID | 2.70 | 1 | ZID |

**4TZK** was selected: highest resolution, single chain, NAD cofactor present, and
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

Centred on the centroid of the co-crystallised inhibitor, [8.79, 32.47, 60.41], with edge
lengths [22.5, 22.0, 22.0] Å (ligand span + 5 Å padding, floored at 22 Å).

The pocket lining was characterised by distance: 22 residues fall within
4.5 Å of ligand 641, of which 20
are apolar or aromatic. Both canonical InhA hydrogen-bonding anchors — **Tyr158** and **Thr196** —
are present, confirming the box covers the substrate-binding site rather than a surface
groove.

## 2. Redocking control

Empirical docking scores are interpretable only if the setup reproduces a known experimental
binding mode. The co-crystallised ligand was extracted, re-protonated, re-prepared as a
flexible ligand, and redocked into the prepared receptor (exhaustiveness 32, 20 poses).

| quantity | value |
|---|---|
| crystal-pose Vina score | -9.38 kcal/mol |
| best redocked score | -9.78 kcal/mol |
| **RMSD, top-ranked pose** | **0.23 Å** |
| best RMSD among 20 poses | 0.23 Å |
| centroid displacement | 0.10 Å |

RMSD note: the docked and crystal poses do not share atom ordering, and a PDBQT-derived
molecule carries no implicit-hydrogen counts, so RDKit's template and substructure-match
routes both fail on this pair. The reported value is an **element-constrained optimal
assignment** (Hungarian algorithm) over heavy atoms in absolute coordinates, with no
realignment — a measurement that cannot flatter the result through superposition.

A top-ranked-pose RMSD of 0.23 Å is well inside the ≤ 2 Å convention for a
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
| ChEMBL records retrieved | 4,450 |
| no data-validity flag | 4,228 |
| exact relation (=) | 3,386 |
| potency endpoint | 2,521 |
| usable units | 2,521 |
| has pChEMBL value | 2,107 |
| pChEMBL in [3,12] | 2,107 |
| unique molecules (aggregated) | 1,417 |

Molecules were standardised with RDKit `rdMolStandardize`: largest-fragment selection
(drops salts and solvates), neutralisation, canonical SMILES. Replicate measurements of the
same molecule against the same target were aggregated to the **median** — whole-cell MICs
are pooled across laboratories and the distribution carries outliers — with the standard
deviation across replicates retained.

**Assay noise floor.** 324 molecules have replicate measurements; their
median inter-assay SD is **0.23 log units**. No model should be expected
to predict below this, and it is the correct scale against which to read the errors in §6.

Final modelling sets: **355** molecules with InhA enzyme potency
(162 active at pChEMBL ≥ 6) and **1097** with whole-cell potency
(283 active). 1,417 unique structures in total; the
difference is molecules measured in both assays.

## 4. Target engagement translates to cell kill

35 molecules were measured in both assays. Their potencies correlate at
**Spearman ρ = 0.91** (p = 4e-14). This is the mechanistic
anchor of the project: InhA inhibition does not merely correlate with an enzymatic readout,
it tracks growth inhibition of the organism.

## 5. QSAR models

**Features.** Morgan r=2, 2048 bits concatenated with 17 physicochemical
descriptors: MW, cLogP, TPSA, HBD, HBA, RotB, RingCount, AromRings, FracCsp3, HeavyAtoms, FormalCharge, QED, BertzCT, NumHeteroatoms, StereoCentres, Lipinski_viol, Veber_ok. The last two are rule-based flags
(Lipinski violation count, Veber pass/fail) rather than continuous properties.

**Labels.** Binary at pChEMBL ≥ 6.0 (1 µM) for both endpoints.

**Splits.** Two, deliberately compared:

- *Scaffold split* — Bemis–Murcko scaffolds; the **largest scaffold groups go to train** and
  rare/singleton scaffolds to test.
  0 scaffolds are shared across the
  enzyme-set boundary and
  0 across the whole-cell boundary
  (152 and
  444 distinct scaffolds respectively).
- *Random split* — same size, uniformly sampled. This leaks whole chemical series across the
  boundary: 29 and
  89 scaffolds appear on both sides for the
  enzyme and whole-cell sets respectively.

The split direction matters and is easy to invert. Assigning the largest scaffold groups to
*test* makes the test set a congeneric series and yields a scaffold-split score **higher**
than the random-split score — a diagnostic that the split is backwards, not that the model
generalises.

**Models.**

- **XGBoost** — 400 trees, depth 6, lr 0.05, subsample 0.8, colsample 0.6, min_child_weight 3, L2 1.5
- **random forest** — 500 trees, min_samples_leaf 2, max_features sqrt, class_weight balanced
- **logistic regression** — L2, C = 0.1, class_weight balanced, on standardised features

**Results (ROC-AUC on held-out test):**

| dataset | model | scaffold | random | inflation |
|---|---|---|---|---|
| InhA enzyme | logistic regression | 0.808 | 0.927 | +0.119 |
| InhA enzyme | random forest | 0.833 | 0.944 | +0.111 |
| InhA enzyme | XGBoost | 0.837 | 0.921 | +0.084 |
| *M. tb* whole-cell | logistic regression | 0.823 | 0.892 | +0.069 |
| *M. tb* whole-cell | random forest | 0.841 | 0.911 | +0.070 |
| *M. tb* whole-cell | XGBoost | 0.793 | 0.868 | +0.076 |

A random split overstates ROC-AUC by
**0.07–0.12** on this data (median
0.08). Every result quoted outside this table is the scaffold-split
number.

## 6. Applicability domain

Test compounds (scaffold split, random forest) were binned by maximum Tanimoto similarity to
any training compound. Bins holding fewer than 10 compounds are flagged: a ROC-AUC computed
on five compounds is not a measurement.

| dataset | Tanimoto bin | n | actives | ROC_AUC | mean abs. err | interpretable |
|---|---|---|---|---|---|---|
| InhA enzyme | 0.00-0.25 | 9 | 4 | 0.450 | 0.509 | **no (n < 10)** |
| InhA enzyme | 0.25-0.35 | 5 | 3 | 1.000 | 0.406 | **no (n < 10)** |
| InhA enzyme | 0.35-0.50 | 13 | 8 | 0.900 | 0.419 | yes |
| InhA enzyme | 0.50-1.01 | 61 | 39 | 0.873 | 0.253 | yes |
| *M. tb* whole-cell | 0.00-0.25 | 40 | 8 | 0.605 | 0.433 | yes |
| *M. tb* whole-cell | 0.25-0.35 | 34 | 4 | 0.475 | 0.426 | yes |
| *M. tb* whole-cell | 0.35-0.50 | 28 | 3 | 0.933 | 0.358 | yes |
| *M. tb* whole-cell | 0.50-1.01 | 172 | 41 | 0.894 | 0.286 | yes |

Accuracy is a function of chemical distance, not a fixed property of the model. Above
Tanimoto 0.35 the interpretable bins score 0.87–0.93; below 0.35 they fall
to 0.47–0.60 — at, or barely above, chance. The trend is not strictly
monotonic (the 0.35–0.50 bins score marginally above the 0.50+ bins, on 13 and 28 compounds),
so the useful reading is the split between the two low-similarity bins and everything else,
not a smooth curve. Predictions on novel chemotypes should not be trusted, and this table is
the reason.

## 7. Docking screen

Library: 647 molecules (347 with InhA enzyme potency,
300 whole-cell), of which 305 are measured actives at
pChEMBL ≥ 6. One ETKDG conformer per molecule, MMFF-minimised, Meeko-prepared.

| dataset | molecules | active (pChEMBL ≥ 6) | % active |
|---|---|---|---|
| InhA enzyme | 347 | 155 | 45 % |
| *M. tb* whole-cell | 300 | 150 | 50 % |

Two selection decisions are worth stating, because both change what the screen can show.

**The whole-cell arm is class-balanced.** Molecules are drawn one-per-scaffold within each
class separately — 150 actives and 150 inactives — rather than
sampled from the dataset as a whole. An earlier version sampled without regard to class and
the arm came out 100 % active, which makes discrimination untestable on that dataset: ROC-AUC
is undefined when one class is empty. Selecting actives by potency instead would be worse
than a bug — it would place the most potent compounds on one side of the comparison and
inflate the apparent separation.

**Molecules outside a drug-like envelope are excluded** (MW ≤ 600 Da,
≤ 12 rotatable bonds; the library's actual maxima are 597 Da
and 12). Vina's scoring function is parameterised on drug-like
complexes, and at fixed exhaustiveness the conformational search under-samples very flexible
ligands, so a score for the most extreme curated molecule (1,876 Da,
62 rotatable bonds) is not comparable to one for a 350 Da fragment.
The envelope removes 71 of the curated molecules. The cut is on ligand
properties only, never on the score, and is applied before docking, so it cannot select on
the outcome. Selection and the single-conformer
limitation are implemented in `src/06_build_library.py`.

## 8. What the docking screen shows

| dataset | n | actives | ROC-AUC | Spearman ρ | EF 5 % | ceiling |
|---|---|---|---|---|---|---|
| InhA enzyme | 342 | 150 | 0.65 | -0.23 | 1.34 | 2.28 |
| *M. tb* whole-cell | 298 | 148 | 0.58 | -0.14 | 0.72 | 2.01 |

Ranking by Vina score separates measured actives from inactives only slightly better than
chance on the enzyme arm and less than that on the whole-cell arm. Rank correlation with
measured potency is weak and negative in the expected direction (a more negative score
should mean a more potent compound), but far too weak to rank a series. Enrichment in the
top 5 % of the ranked list is below the attainable ceiling on both arms, and on the
whole-cell arm it is below 1.0 - selecting the top-scoring 5 % of that library would return
fewer actives than picking at random.

**The comparison that matters is on the same molecules.** Docking and the QSAR model were
scored on the held-out scaffold-split test compounds, which neither approach was fitted to:

| dataset | held-out n | docking ROC-AUC | QSAR ROC-AUC |
|---|---|---|---|
| InhA enzyme | 82 | 0.59 | 0.83 |
| *M. tb* whole-cell | 140 | 0.63 | 0.84 |

On identical molecules the ligand-based model separates actives from inactives far better
than the physics-based score. This is the expected result for a target with this much
measured data - a model that has seen 2,107 curated activity records for
1,417 molecules is using information that a single rigid-receptor docking run
does not have. The screen is reported because a negative result on a well-posed comparison
is a result: it bounds what the docking component can contribute, and it is the reason the
prediction endpoint is the QSAR model rather than the docking score.

Three specific limitations bound the docking number, and none of them are ruled out here.
One conformer per ligand under-samples flexible molecules. The receptor is a single rigid
crystal structure, so induced fit is not modelled. And Vina's score correlates with
molecular size (Spearman rho vs MW is reported in figure 3f), which inflates the apparent
performance of any list ranked by raw score. A size-corrected score, ensemble docking across
the 10 InhA structures screened, or rescoring the retained poses would each test
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
