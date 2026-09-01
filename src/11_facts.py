"""Step 11 - derive every headline number from the saved tables and inject them into the docs.

The README and METHODS.md must never carry a hand-typed number: an edit to a filter or a
hyperparameter would silently leave the prose claiming an old result. This script is the
single source of truth. It writes docs/facts.json and rewrites the block between the
FACTS markers in README.md.

Run it after 04 (QSAR) and, once docking has finished, after 09 (docking analysis).
"""
import json
import os, os, re
import numpy as np, pandas as pd
from scipy import stats

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
START, END = "<!-- FACTS:START -->", "<!-- FACTS:END -->"
AD_START, AD_END = "<!-- AD:START -->", "<!-- AD:END -->"
FIG_START, FIG_END = "<!-- FIGURES:START -->", "<!-- FIGURES:END -->"


def facts():
    R = pd.read_csv(f"{D}/data/metrics.csv")
    AG = pd.read_csv(f"{D}/data/activity_aggregated.csv")
    AD = pd.read_csv(f"{D}/data/applicability_domain.csv")
    SP = json.load(open(f"{D}/data/splits.json"))
    CAS = json.load(open(f"{D}/data/curation_cascade.json"))["cascade"]
    V = json.load(open(f"{D}/dock/redock_validation.json"))
    BOX = json.load(open(f"{D}/dock/box.json"))
    LIB = pd.read_csv(f"{D}/dock/library.csv")
    FT = pd.read_csv(f"{D}/data/features.csv")
    SM = json.load(open(f"{D}/data/structure_meta.json"))

    f = {}
    f["cascade"] = [[k, int(v)] for k, v in CAS]
    f["records_curated"] = int([v for k, v in CAS if "pChEMBL in" in k][0])
    f["molecules_unique"] = int(AG.std_smiles.nunique())
    f["rows_mol_by_dataset"] = int(len(AG))
    for ds, tag in [("InhA_enzyme", "inha"), ("Mtb_whole_cell", "wc")]:
        g = AG[AG.dataset == ds]
        f[f"n_{tag}"] = int(len(g))
        f[f"actives_{tag}"] = int((g.pchembl >= 6).sum())
    rep = AG[AG.n_meas > 1]
    f["n_replicated"] = int(len(rep))
    f["noise_floor_sd"] = round(float(rep.pchembl_sd.median()), 3)

    # enzyme -> whole-cell translation
    p = AG.pivot_table(index="std_smiles", columns="dataset", values="pchembl").dropna()
    rho, pv = stats.spearmanr(p.InhA_enzyme, p.Mtb_whole_cell)
    f["translation"] = dict(spearman_rho=round(float(rho), 3), n=int(len(p)), p_value=float(pv))

    # QSAR: best model per dataset on the scaffold split, with its random-split counterpart
    f["qsar"] = {}
    for ds in R.dataset.unique():
        sc = R[(R.dataset == ds) & (R.split == "scaf")].sort_values("ROC_AUC").iloc[-1]
        rd = R[(R.dataset == ds) & (R.split == "rand") & (R.model == sc.model)].iloc[0]
        f["qsar"][ds] = dict(
            best_model=sc.model,
            scaffold_roc_auc=round(float(sc.ROC_AUC), 3), random_roc_auc=round(float(rd.ROC_AUC), 3),
            scaffold_pr_auc=round(float(sc.PR_AUC), 3), scaffold_mcc=round(float(sc.MCC), 3),
            n_train=int(sc.n_train), n_test=int(sc.n_test), test_actives=int(sc.test_actives),
            ef10=round(float(sc.EF10), 2), ef_max=round(float(sc.EF_max), 2),
            n_scaffolds=int(SP["splits"][ds]["n_scaffolds"]),
            leaked_scaffolds_random=int(SP["splits"][ds]["leaked_scaffolds_random"]))
    piv = R.pivot_table(index=["dataset", "model"], columns="split", values="ROC_AUC")
    infl = (piv["rand"] - piv["scaf"])
    f["inflation"] = dict(min=round(float(infl.min()), 3), max=round(float(infl.max()), 3),
                          median=round(float(infl.median()), 3))

    # applicability domain: near vs far, ignoring bins too small to interpret
    ok = AD[AD.n >= 10]
    near = ok[~ok.tanimoto_bin.str.startswith(("0.00", "0.25"))]   # >= 0.35 similarity
    far = ok[ok.tanimoto_bin.str.startswith(("0.00", "0.25"))]      # < 0.35 similarity
    f["applicability_domain"] = dict(
        near_roc_auc=[round(float(v), 3) for v in near.ROC_AUC],
        far_roc_auc=[round(float(v), 3) for v in far.ROC_AUC],
        n_bins_too_small=int((AD.n < 10).sum()),
        table=AD.round(3).to_dict("records"))

    # docking setup
    f["docking"] = dict(
        receptor=SM["receptor"] if "receptor" in SM else "4TZK",
        rmsd_top=V["rmsd_top"], rmsd_best=V["rmsd_best"],
        crystal_score=V["crystal_score"], redock_top=V["redock_top5"][0],
        box_size=BOX["size"], box_center=[round(c, 2) for c in BOX["center"]],
        library_n=int(len(LIB)),
        library_inha=int((LIB.dataset == "InhA_enzyme").sum()),
        library_wc=int((LIB.dataset == "Mtb_whole_cell").sum()),
        library_actives=int((LIB.pchembl >= 6).sum()),
        library_mw_max=float(round(LIB.MW.max(), 1)),
        library_rotb_max=int(LIB.RotB.max()),
        # the most extreme curated molecule the envelope excludes - quoted in METHODS as the
        # concrete reason for the cut, so it has to track the data
        excluded_max_mw=float(round(FT.MW.max(), 0)),
        excluded_max_rotb=int(FT.RotB.max()),
        n_excluded=int(((FT.MW > 600) | (FT.RotB > 12)).sum()),
        # per-dataset balance: a one-class arm makes discrimination untestable there
        library_balance={ds: dict(n=int(len(g)), actives=int((g.pchembl >= 6).sum()))
                         for ds, g in LIB.groupby("dataset")})
    sc = f"{D}/data/docking_joined.csv"
    if os.path.exists(sc):
        S = pd.read_csv(sc)
        f["docking"]["n_scored"] = int(len(S))
        if {"vina_kcal", "pchembl"} <= set(S.columns):
            r, pvv = stats.spearmanr(S.vina_kcal, S.pchembl)
            f["docking"]["score_vs_potency_rho"] = round(float(r), 3)
            f["docking"]["score_vs_potency_p"] = float(pvv)
    pf = f"{D}/data/docking_performance.json"
    if os.path.exists(pf):
        P = json.load(open(pf))
        for ds, d in P["by_dataset"].items():
            f["docking"][f"{ds}_roc_auc"] = round(float(d["roc_auc"]), 3)
            f["docking"][f"{ds}_ef5"] = round(float(d["ef5"]), 2)
            f["docking"][f"{ds}_ef_max"] = round(float(d["ef_max"]), 2)
            f["docking"][f"{ds}_rho"] = round(float(d["spearman"]), 3)
        for ds, d in P["held_out_head_to_head"].items():
            f["docking"][f"{ds}_heldout_n"] = int(d["docking"]["n"])
            f["docking"][f"{ds}_heldout_dock_auc"] = round(float(d["docking"]["roc_auc"]), 3)
            f["docking"][f"{ds}_heldout_qsar_auc"] = round(float(d["qsar"]["roc_auc"]), 3)
    lk = f"{D}/data/leakage_decomposition.csv"
    if os.path.exists(lk):
        L = pd.read_csv(lk)
        g = L.dropna(subset=["delta_leakage"])
        rf = L[L.model == "RandomForest"]
        f["leakage"] = {
            "n_combinations": int(len(g)),
            "n_positive": int((g.delta_leakage > 0).sum()),
            "median_delta": round(float(g.delta_leakage.median()), 3),
            "median_residual": round(float(g.delta_residual.median()), 3),
        }
        for _, r in rf.iterrows():
            f["leakage"][f"{r.dataset}_{r.group}_auc"] = round(float(r.roc_auc), 2)
            f["leakage"][f"{r.dataset}_{r.group}_n"] = int(r.n)

    f["descriptors_n"] = len(SP["descriptor_cols"])
    f["fingerprint"] = f"Morgan r={SP['radius']}, {SP['nbits']} bits"
    return f


def headline_md(f):
    q, dk = f["qsar"], f["docking"]
    ie, wc = q["InhA_enzyme"], q["Mtb_whole_cell"]
    rows = [
        ("Redocking control (top-ranked pose vs. crystal)",
         f"**{f['docking']['rmsd_top']:.2f} Å** RMSD"),
        ("Compounds in the docking library", f"{f['docking']['library_n']}"),
        ("Curated bioactivity records",
         f"{f['records_curated']:,} → {f['molecules_unique']:,} unique molecules"),
        ("Enzyme → whole-cell potency correlation",
         f"**Spearman ρ = {f['translation']['spearman_rho']:.2f}** "
         f"(n = {f['translation']['n']})"),
        (f"InhA classifier, scaffold split ({ie['best_model']})",
         f"ROC-AUC **{ie['scaffold_roc_auc']:.2f}** "
         f"(random split: {ie['random_roc_auc']:.2f})"),
        (f"Whole-cell classifier, scaffold split ({wc['best_model']})",
         f"ROC-AUC **{wc['scaffold_roc_auc']:.2f}** "
         f"(random split: {wc['random_roc_auc']:.2f})"),
        ("Experimental noise floor (median SD across replicate assays)",
         f"{f['noise_floor_sd']:.2f} log units"),
        # the table claimed a docking library without ever saying what the screen found
        ("Docking vs. QSAR on the same held-out molecules",
         f"ROC-AUC **{dk['InhA_enzyme_heldout_dock_auc']:.2f}** / "
         f"**{dk['Mtb_whole_cell_heldout_dock_auc']:.2f}** "
         f"vs. **{dk['InhA_enzyme_heldout_qsar_auc']:.2f}** / "
         f"**{dk['Mtb_whole_cell_heldout_qsar_auc']:.2f}**"),
    ]
    md = ["| | value |", "|---|---|"] + [f"| {k} | {v} |" for k, v in rows]
    md += ["", f"A random train/test split leaks chemical series across the boundary and "
               f"inflates ROC-AUC by **{f['inflation']['min']:.2f}–{f['inflation']['max']:.2f}** "
               f"(median {f['inflation']['median']:.2f}) across all six model × dataset "
               f"combinations. Every number quoted as a result here is the scaffold-split "
               f"number."]
    return "\n".join(md)


def figures_md(f):
    """Embed the figure set. Captions carry the number, not just the title, so the
    landing page states a result even if the reader never opens the figure."""
    dk, q = f["docking"], f["qsar"]
    ie, wc = q["InhA_enzyme"], q["Mtb_whole_cell"]
    figs = [
        ("figures/fig1_data_and_target.png", "Dataset and target",
         f"Curated bioactivity landscape ({f['records_curated']:,} records → "
         f"{f['molecules_unique']:,} unique molecules), enzyme→whole-cell potency "
         f"translation (Spearman ρ = {f['translation']['spearman_rho']:.2f}), and the "
         f"InhA binding site used for docking."),
        ("figures/fig2_qsar_validation.png", "QSAR validation",
         f"Scaffold-split performance ({ie['best_model']} {ie['scaffold_roc_auc']:.2f} on "
         f"InhA, {wc['best_model']} {wc['scaffold_roc_auc']:.2f} whole-cell), the "
         f"random-split inflation gap, and accuracy as a function of chemical distance."),
        ("figures/fig3_docking.png", "Docking screen",
         f"Score distributions, score vs. measured potency, ROC, and the head-to-head on "
         f"held-out molecules: docking "
         f"{dk['InhA_enzyme_heldout_dock_auc']:.2f}/"
         f"{dk['Mtb_whole_cell_heldout_dock_auc']:.2f} vs. QSAR "
         f"{dk['InhA_enzyme_heldout_qsar_auc']:.2f}/"
         f"{dk['Mtb_whole_cell_heldout_qsar_auc']:.2f} ROC-AUC."),
    ]
    md = []
    for path, title, cap in figs:
        md += [f"### {title}", "", f"![{title}]({path})", "", f"*{cap}*", ""]
    return "\n".join(md).rstrip()


def ad_md(f):
    ad = f["applicability_domain"]
    near, far = ad["near_roc_auc"], ad["far_roc_auc"]
    return (f"- **Accuracy is a function of chemical distance.** Test compounds at Tanimoto "
            f"similarity >= 0.35 to the training set score ROC-AUC "
            f"{min(near):.2f}–{max(near):.2f}; below 0.35 the same model falls to "
            f"{min(far):.2f}–{max(far):.2f} — at, or barely above, chance. Predictions on\n"
            f"  genuinely novel chemotypes are not trustworthy, and Figure 2d shows it rather\n"
            f"  than hiding it. ({ad['n_bins_too_small']} similarity bins hold fewer than 10\n"
            f"  compounds; they are hatched in the figure and excluded from these ranges "
            f"rather than\n  quoted as results.)")


def inject(path, start, end, body, label):
    s = open(path).read()
    if start in s and end in s:
        pre, rest = s.split(start, 1)
        _, post = rest.split(end, 1)
        open(path, "w").write(f"{pre}{start}\n{body}\n{end}{post}")
        print(f"{label} updated")
    else:
        print(f"WARNING: {start} / {end} markers not found in {path}")


def main():
    f = facts()
    os.makedirs(f"{D}/docs", exist_ok=True)
    json.dump(f, open(f"{D}/docs/facts.json", "w"), indent=1)

    rd = f"{D}/README.md"
    inject(rd, START, END, headline_md(f), "README headline table")
    inject(rd, AD_START, AD_END, ad_md(f), "README applicability-domain bullet")
    inject(rd, FIG_START, FIG_END, figures_md(f), "README figure gallery")
    print(json.dumps({k: v for k, v in f.items()
                      if k not in ("cascade", "applicability_domain", "qsar")}, indent=1))


if __name__ == "__main__":
    main()
