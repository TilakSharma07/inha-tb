"""Step 9 - regenerate the publication figures from saved tables only.

Runs without re-doing the docking screen: everything is read from data/*.csv,
data/*.json and dock/*.json. Two figures are produced:

  fig1_data_and_target.png  - curation cascade, potency distributions, assay noise floor,
                              structure selection, binding-pocket contacts, chemical space
  fig2_qsar_validation.png  - ROC, split comparison, precision-recall, applicability
                              domain, enzyme-vs-cell correlation, chemotype separation
"""
import json, os
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FIG = f"{D}/figures"
FOCAL, COMP, LIGHT, GREY = "#1f4e79", "#c1660a", "#4a7ba7", "#8c8c8c"
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "PHE", "MET", "TRP", "PRO", "GLY", "TYR"}

mpl.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "font.size": 7,
    "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.2,
    "ytick.labelsize": 6.2, "legend.fontsize": 6, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 1.3,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans", "figure.constrained_layout.use": False,
})


def panel_letter(ax, L, dx=-0.18, dy=1.02):
    ax.text(dx, dy, L, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


def figure1():
    CAS = json.load(open(f"{D}/data/curation_cascade.json"))["cascade"]
    # drop stages that removed no records - an identical-length bar carries no information
    CAS = [kv for i, kv in enumerate(CAS)
           if i == 0 or i == len(CAS) - 1 or kv[1] != CAS[i - 1][1]]
    AG = pd.read_csv(f"{D}/data/activity_aggregated.csv")
    FT = pd.read_csv(f"{D}/data/features.csv")
    SM = json.load(open(f"{D}/data/structure_meta.json"))
    POK = json.load(open(f"{D}/dock/pocket_residues.json"))
    lib = pd.read_csv(f"{D}/dock/library.csv")

    fig = plt.figure(figsize=(7.2, 6.4))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.52, wspace=0.42)
    a, b, c = [fig.add_subplot(gs[0, i]) for i in range(3)]
    d, e, f = [fig.add_subplot(gs[1, i]) for i in range(3)]

    # a - curation cascade
    yv, vals = np.arange(len(CAS))[::-1], [v for _, v in CAS]
    a.barh(yv, vals, color=[FOCAL] * (len(CAS) - 1) + [COMP], height=0.62)
    for y, v in zip(yv, vals):
        a.text(v + 90, y, f"{v:,}", va="center", fontsize=5.8, color="#333")
    a.set_yticks(yv)
    a.set_yticklabels([k.replace("ChEMBL records retrieved", "ChEMBL records\nretrieved")
                        .replace(" (=)", "\n(=)")
                        .replace("no data-validity flag", "no data-validity\nflag")
                        .replace("unique molecules\n(median-aggregated)",
                                 "unique molecules\n(aggregated)")
                       for k, _ in CAS], fontsize=5.5)
    a.set_xlim(0, 6000); a.set_xlabel("records", fontsize=6.0)
    a.set_title("Curation cascade", loc="left")

    # b - potency distributions
    for ds, col, lab in [("InhA_enzyme", FOCAL, "InhA enzyme  p(IC$_{50}$/K$_i$)"),
                         ("Mtb_whole_cell", LIGHT, "M. tb cells  p(MIC)")]:
        v = AG[AG.dataset == ds].pchembl
        b.hist(v, bins=np.arange(3.25, 9.75, 0.4), histtype="stepfilled", alpha=0.62,
               color=col, label=f"{lab}\nn={len(v)}", lw=0.8, edgecolor=col)
    b.axvline(6.0, color=GREY, ls="--", lw=0.8)
    b.text(5.92, b.get_ylim()[1] * 0.42, "active\nthreshold", fontsize=5.4, color=GREY,
           va="center", ha="right")
    b.set_xlabel("potency (pChEMBL)"); b.set_ylabel("molecules")
    b.set_title("Potency distributions", loc="left")
    b.legend(frameon=False, fontsize=5.2, loc="upper right", handlelength=1.0,
             labelspacing=0.5, bbox_to_anchor=(1.02, 1.02), borderaxespad=0.0)
    b.set_ylim(top=b.get_ylim()[1] * 1.42)

    # c - replicate spread = experimental noise floor
    sd = AG[AG.n_meas > 1].pchembl_sd.dropna()
    c.hist(sd, bins=np.arange(0, 2.6, 0.12), color=FOCAL, alpha=0.8)
    med = sd.median()
    c.axvline(med, color=COMP, lw=1.2)
    c.text(med + 0.08, c.get_ylim()[1] * 0.88, f"median\n{med:.2f} log units",
           fontsize=5.6, color=COMP)
    c.set_xlabel("SD across replicate assays (log units)"); c.set_ylabel("molecules")
    c.set_title("Experimental noise floor", loc="left")

    # d - structure selection by resolution
    S = pd.DataFrame([dict(pdb=k, res=v.get("resolution_A")) for k, v in SM.items()])
    S = S[S.res.notna()].sort_values("res")
    yv = np.arange(len(S))[::-1]
    d.barh(yv, S.res, color=[COMP if p == "4TZK" else FOCAL for p in S.pdb], height=0.6)
    for y, (_, r) in zip(yv, S.iterrows()):
        tag = "  \u2190 receptor" if r.pdb == "4TZK" else ""
        d.text(r.res + 0.05, y, f"{r.res:.2f}{tag}", va="center", fontsize=5.5,
               color=COMP if r.pdb == "4TZK" else "#333")
    d.set_yticks(yv); d.set_yticklabels(S.pdb, fontsize=5.8)
    d.set_xlabel("resolution (\u00c5)  \u2014 lower is better"); d.set_xlim(0, 3.6)
    d.set_title("InhA structures screened", loc="left")

    # e - binding-pocket contacts
    pk = pd.DataFrame(POK).sort_values("min_dist")
    yv = np.arange(len(pk))[::-1]
    e.barh(yv, pk.min_dist, height=0.68,
           color=[FOCAL if r.resn in HYDROPHOBIC else COMP for _, r in pk.iterrows()])
    e.set_yticks(yv)
    e.set_yticklabels([f"{r.resn.title()}{r.resi}" for _, r in pk.iterrows()], fontsize=4.6)
    e.set_xlabel("min. distance to inhibitor 641 (\u00c5)"); e.set_xlim(2.5, 4.7)
    e.set_title("Binding-pocket contacts", loc="left")
    nh = int(pk.resn.isin(HYDROPHOBIC).sum())
    e.text(1.0, 1.005, f"{nh}/{len(pk)} residues hydrophobic", transform=e.transAxes,
           ha="right", va="bottom", fontsize=5.0, color=GREY)

    # f - chemical space: docked subset vs curated set
    f.scatter(FT.MW, FT.cLogP, s=3.5, c="#c9c9c9", alpha=0.55, lw=0,
              label=f"curated set ({len(FT)})")
    f.scatter(lib.MW, lib.cLogP, s=4.5, c=COMP, alpha=0.75, lw=0, label=f"docked ({len(lib)})")
    f.axvline(500, color=GREY, ls=":", lw=0.7); f.axhline(5, color=GREY, ls=":", lw=0.7)
    f.set_xlabel("molecular weight (Da)"); f.set_ylabel("cLogP")
    f.set_title("Docked subset vs. curated set", loc="left")
    f.legend(frameon=False, fontsize=5.2, loc="upper left")

    for ax, L in zip([a, b, c, d, e, f], "abcdef"):
        panel_letter(ax, L, dx=-0.30 if ax is b else -0.18)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/fig1_data_and_target.png", dpi=300, bbox_inches="tight")
    print("wrote fig1_data_and_target.png")
    return fig


def figure2():
    R = pd.read_csv(f"{D}/data/metrics.csv")
    CU = json.load(open(f"{D}/data/curves.json"))
    AD = pd.read_csv(f"{D}/data/applicability_domain.csv")
    AG = pd.read_csv(f"{D}/data/activity_aggregated.csv")
    FT = pd.read_csv(f"{D}/data/features.csv")
    SP = json.load(open(f"{D}/data/splits.json"))
    X = np.load(f"{D}/data/fp_matrix.npy")

    fig = plt.figure(figsize=(7.2, 6.4))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.52, wspace=0.42)
    a, b, c = [fig.add_subplot(gs[0, i]) for i in range(3)]
    d, e, f = [fig.add_subplot(gs[1, i]) for i in range(3)]

    # a - ROC curves, both datasets x both splits
    for key, col, ls, lab in [("InhA_enzyme|scaf", FOCAL, "-", "InhA, scaffold"),
                              ("InhA_enzyme|rand", COMP, "-", "InhA, random"),
                              ("Mtb_whole_cell|scaf", FOCAL, "--", "whole-cell, scaffold"),
                              ("Mtb_whole_cell|rand", COMP, "--", "whole-cell, random")]:
        cv = CU[key]
        ds, sp = key.split("|")
        auc = R[(R.dataset == ds) & (R.split == sp) & (R.model == "XGBoost")].ROC_AUC.iloc[0]
        a.plot(cv["fpr"], cv["tpr"], color=col, ls=ls, label=f"{lab} ({auc:.2f})")
    a.plot([0, 1], [0, 1], ":", color=GREY, lw=0.8)
    a.set_xlim(-0.02, 1.02); a.set_ylim(-0.02, 1.02)
    a.set_xlabel("False positive rate"); a.set_ylabel("True positive rate")
    a.set_title("Random splits overstate accuracy", loc="left")
    a.legend(frameon=False, fontsize=5.4, loc="lower right")

    # b - same model, two splits
    MS = {"XGBoost": "XGBoost", "RandomForest": "random forest", "LogReg(baseline)": "logistic reg."}
    MORDER = ["LogReg(baseline)", "RandomForest", "XGBoost"]
    order = [("InhA_enzyme", m) for m in MORDER] + [("Mtb_whole_cell", m) for m in MORDER]
    piv = R.pivot_table(index=["dataset", "model"], columns="split", values="ROC_AUC")
    yp = list(range(len(order), 0, -1))
    for (ds, mo), y in zip(order, yp):
        row = piv.loc[(ds, mo)]
        b.plot([row["scaf"], row["rand"]], [y, y], "-", color=GREY, lw=1.0, zorder=1)
        b.scatter(row["scaf"], y, s=42, color=FOCAL, zorder=3)
        b.scatter(row["rand"], y, s=42, color=COMP, marker="D", zorder=3)
    b.set_yticks(yp); b.set_yticklabels([MS[m] for _, m in order], fontsize=6.0)
    b.set_ylim(-1.0, len(order) + 0.6); b.set_xlim(0.62, 1.0)
    b.set_xlabel("ROC-AUC on held-out test set")
    b.set_title("Same model, two splits", loc="left")
    b.axhline(3.5, color=GREY, lw=0.6, alpha=0.6)
    b.text(0.03, np.mean(yp[:3]), "InhA\nenzyme", transform=b.get_yaxis_transform(),
           fontsize=6.0, ha="left", va="center", style="italic", color=FOCAL)
    b.text(0.03, np.mean(yp[3:]), "M. tb\ncells", transform=b.get_yaxis_transform(),
           fontsize=6.0, ha="left", va="center", color=LIGHT)
    b.scatter([], [], s=42, color=FOCAL, label="scaffold split (honest)")
    b.scatter([], [], s=42, color=COMP, marker="D", label="random split (leaky)")
    b.legend(frameon=False, fontsize=5.8, loc="lower left", borderaxespad=0.1,
             handletextpad=0.4, labelspacing=0.22)

    # c - precision-recall, scaffold split
    for key, col, lab in [("InhA_enzyme|scaf", FOCAL, "InhA enzyme"),
                          ("Mtb_whole_cell|scaf", LIGHT, "whole-cell")]:
        cv = CU[key]
        ds = key.split("|")[0]
        ap = R[(R.dataset == ds) & (R.split == "scaf") & (R.model == "XGBoost")].PR_AUC.iloc[0]
        c.plot(cv["rec"], cv["prec"], color=col, label=f"{lab} (AP {ap:.2f})")
        c.axhline(cv["base"], color=col, ls=":", lw=0.8)
    c.set_xlim(-0.02, 1.02); c.set_ylim(-0.02, 1.05)
    c.set_xlabel("Recall"); c.set_ylabel("Precision")
    c.set_title("Precision\u2013recall, scaffold split", loc="left")
    c.legend(frameon=False, fontsize=5.4, loc="upper right")
    c.text(0.03, 0.06, "dotted = random baseline\n(prevalence)", transform=c.transAxes,
           fontsize=5.0, color=GREY)

    # d - applicability domain
    bins = AD.tanimoto_bin.unique().tolist()
    w, xx = 0.38, np.arange(len(bins))
    for i, (ds, col, lab) in enumerate([("InhA_enzyme", FOCAL, "InhA enzyme"),
                                        ("Mtb_whole_cell", LIGHT, "whole-cell")]):
        sub = AD[AD.dataset == ds].set_index("tanimoto_bin").reindex(bins)
        small = (sub.n < 10).fillna(False).to_numpy()
        d.bar(xx + (i - 0.5) * w, sub.ROC_AUC, w, color=col, label=lab)
        d.bar(xx[small] + (i - 0.5) * w, sub.ROC_AUC[small], w, color="white", alpha=0.55,
              hatch="////", edgecolor=col, lw=0.5, zorder=2.5)
        for x, (auc, n) in zip(xx + (i - 0.5) * w, zip(sub.ROC_AUC, sub.n)):
            if auc == auc:
                d.text(x, 0.315, f"n={int(n)}", ha="center", va="bottom", rotation=90,
                       fontsize=4.6, color="white")
    d.axhline(0.5, color=GREY, ls=":", lw=0.9)
    d.text(len(bins) - 0.55, 0.505, "random", fontsize=5.0, color=GREY, va="bottom", ha="right")
    d.set_xticks(xx); d.set_xticklabels([("<0.25" if b.startswith("0.00") else
                        "\u22650.50" if b.endswith("1.01") else
                        b.replace("-", "\u2013\n")) for b in bins], fontsize=5.6)
    d.set_ylim(0.3, 1.06); d.set_ylabel("ROC-AUC")
    d.set_xlabel("Max Tanimoto similarity to training set")
    d.set_title("Accuracy vs chemical distance", loc="left")
    d.text(0.985, 0.88, "hatched: n < 10", transform=d.transAxes,
           ha="right", va="top", fontsize=5.0, color=GREY)
    d.legend(frameon=False, fontsize=5.2, loc="upper left", handlelength=1.0,
             borderaxespad=0.2, labelspacing=0.25)

    # e - enzyme potency vs whole-cell potency
    p = AG.pivot_table(index="std_smiles", columns="dataset", values="pchembl").dropna()
    rho, pv = stats.spearmanr(p.InhA_enzyme, p.Mtb_whole_cell)
    lo = min(p.InhA_enzyme.min(), p.Mtb_whole_cell.min()) - 0.5
    hi = max(p.InhA_enzyme.max(), p.Mtb_whole_cell.max()) + 0.5
    e.plot([lo, hi], [lo, hi], ":", color=GREY, lw=0.8)
    e.scatter(p.InhA_enzyme, p.Mtb_whole_cell, s=13, color=LIGHT, alpha=0.85, lw=0)
    e.set_xlim(lo, hi); e.set_ylim(lo, hi)
    e.set_xlabel("InhA enzyme  p(IC$_{50}$/K$_i$)")
    e.set_ylabel("M. tuberculosis cells  p(MIC)")
    e.set_title("Enzyme potency tracks cell kill", loc="left")
    e.text(0.97, 0.20, f"Spearman \u03c1 = {rho:.2f}\nn = {len(p)}\np = {pv:.0e}",
           transform=e.transAxes, ha="right", fontsize=5.8)
    e.text(0.97, 0.05, "dotted = identity", transform=e.transAxes, ha="right",
           fontsize=5.0, color=GREY)

    # f - chemotype separation of the scaffold split (fingerprint PCA)
    sp = SP["splits"]["Mtb_whole_cell"]
    gi = np.array(sp["global_idx"])
    Xi = X[gi].astype(float)
    Xc = Xi - Xi.mean(0)
    U, S_, Vt = np.linalg.svd(Xc, full_matrices=False)
    PC = U[:, :2] * S_[:2]
    tr, te = np.array(sp["scaf_train"]), np.array(sp["scaf_test"])
    f.scatter(PC[tr, 0], PC[tr, 1], s=5, c="#c9c9c9", alpha=0.7, lw=0, label=f"train ({len(tr)})")
    f.scatter(PC[te, 0], PC[te, 1], s=5, c=COMP, alpha=0.8, lw=0, label=f"test ({len(te)})")
    f.set_xlabel("Fingerprint PC1", fontsize=6.0); f.set_ylabel("PC2", fontsize=6.0)
    f.set_xticks([]); f.set_yticks([])
    f.set_title("Test chemotypes sit apart", loc="left")
    f.text(0.03, 0.03, "whole-cell scaffold split", transform=f.transAxes,
           fontsize=5.0, color=GREY)
    f.legend(frameon=False, fontsize=5.4, loc="upper left", handletextpad=0.3,
             borderaxespad=0.2)

    for ax, L in zip([a, b, c, d, e, f], "abcdef"):
        panel_letter(ax, L, dx=-0.26 if ax in (c, d) else -0.18)
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/fig2_qsar_validation.png", dpi=300, bbox_inches="tight")
    print("wrote fig2_qsar_validation.png")
    return fig


if __name__ == "__main__":
    for f in (figure1(), figure2()):
        plt.close(f)
