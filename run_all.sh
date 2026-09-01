#!/usr/bin/env bash
# Full pipeline: ChEMBL query -> curated data -> QSAR -> docking -> figures -> docs.
# Steps 1 and 8 are the slow ones (network pull; ~500 docking runs on CPU).
# Step 1 is optional: the repo ships data/chembl_raw.json.gz, so step 2 runs from a clone.
# Step 8 is resumable - re-running it recovers finished ligands from dock/poses/.
set -euo pipefail
cd "$(dirname "$0")"

python src/01_fetch_chembl.py        # ChEMBL -> data/chembl_raw.json
python src/02_curate.py              # filter cascade + standardise + aggregate replicates
python src/03_features.py            # fingerprints, descriptors, scaffold & random splits
python src/04_train_qsar.py          # QSAR models on both splits + applicability domain
python src/05_prepare_receptor.py    # 4TZK -> receptor.pdbqt, docking box
python src/06_build_library.py       # docking library + 3D conformers
python src/07_redock_validate.py     # CONTROL: redock the crystal ligand, assert RMSD < 2 A
python src/08_dock_library.py        # dock the library (parallel; slow)
python src/09_analyse_docking.py     # docking score vs measured potency
python src/10_figures.py             # regenerate figures 1-2 from saved tables
python src/11_facts.py               # docs/facts.json + README headline table
python src/12_methods_doc.py         # docs/METHODS.md, generated from the tables
python src/13_figure_docking.py      # figure 3, from the docking join
python src/17_hinglish_pdf.py         # plain-language explainer PDF (Hinglish)
python src/15_check_consistency.py    # assert docs and tables state the same numbers
python src/16_check_figures.py       # assert no text collisions / clipped labels
python src/14_package.py             # distributable archive
