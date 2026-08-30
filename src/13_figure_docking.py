"""Step 13 - Figure 3: what the docking screen actually delivered.

Regenerated from data/docking_joined.csv and data/docking_performance.json alone, so it
never depends on live kernel state. Panels:

  a  score distribution, actives vs inactives (the honest first look: heavy overlap)
  b  Vina score vs measured potency, per dataset, with Spearman rho
  c  ROC curves: docking score as a classifier of measured activity
  d  enrichment at the top of the ranked list, against the attainable ceiling
  e  head-to-head on molecules the QSAR model never trained on
  f  score vs molecular weight - the size bias every Vina screen has and few report
"""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FOCAL, COMP, LIGHT, GREY = "#1f4e79", "#c55a11", "#7ba7d0", "#808080"
DS_LABEL = {"InhA_enzyme": "InhA enzyme", "Mtb_whole_cell": "M. tb whole-cell"}


def style():
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300, "font.size": 7,
        "axes.titlesize": 8, "axes.labelsize": 7, "axes.spines.top": False,
        "axes.spines.right": False, "xtick.labelsize": 6.2, "ytick.labelsize": 6.2,
        "legend.fontsize": 6, "axes.linewidth": 0.6, "lines.linewidth": 1.1,
        "axes.titlelocation": "left", "axes.titlepad": 4,
    })


def panel_letter(ax, L):
    ax.text(-0.26, 1.16, L, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="top", ha="left")


def main():
    J = f"{D}/data/docking_joined.csv"
    if not os.path.exists(J):
        sys.exit("data/docking_joined.csv absent - run src/09_analyse_docking.py first")
    df = pd.read_csv(J)
    P = json.load(open(f"{D}/data/docking_performance.json"))
    style()
    fig, AX = plt.subplots(2, 3, figsize=(7.4, 5.0))
    (a, b, c), (d, e, f) = AX
    fig.subplots_adjust(left=0.088, right=0.985, top=0.89, bottom=0.09,
                        wspace=0.46, hspace=0.58)

    # a - score distributions
    bins = np.linspace(df.vina_kcal.min(), df.vina_kcal.max(), 26)
    for lab, k, col in [("inactive", df.active == 0, LIGHT), ("active", df.active == 1, FOCAL)]:
        a.hist(df.vina_kcal[k], bins=bins, alpha=0.72, color=col,
               label=f"{lab} (n={int(k.sum())})", edgecolor="white", linewidth=0.3)
    a.set_xlabel("Vina score (kcal/mol)")
    a.set_ylabel("molecules")
    a.set_title("Actives and inactives overlap")
    a.legend(frameon=False, loc="upper left")

    # b - score vs potency
    for i, (ds, sub) in enumerate(df.groupby("dataset")):
        col = FOCAL if i == 0 else COMP
        b.scatter(sub.vina_kcal, sub.pchembl, s=7, alpha=0.5, color=col,
                  edgecolors="none", label=DS_LABEL.get(ds, ds))
    b.axhline(6.0, color=GREY, ls=":", lw=0.8)
    b.annotate("active threshold", xy=(0.99, 6.0), xycoords=("axes fraction", "data"),
               xytext=(0, 3), textcoords="offset points", ha="right", va="bottom",
               fontsize=5.4, color=GREY)
    b.set_xlabel("Vina score (kcal/mol)")
    b.set_ylabel("measured potency (pChEMBL)")
    b.set_title("Score barely tracks potency")
    lg = b.legend(frameon=False, loc="lower left", scatterpoints=1, handletextpad=0.4)
    for h in lg.legend_handles:
        h.set_alpha(1.0)
        h.set_sizes([22])

    # c - ROC of the docking score
    for i, (ds, sub) in enumerate(df.groupby("dataset")):
        y = sub.active.to_numpy()
        if y.sum() < 5 or y.sum() == len(y):
            continue
        fpr, tpr, _ = roc_curve(y, -sub.vina_kcal)
        auc = P["by_dataset"].get(ds, {}).get("roc_auc", np.nan)
        c.plot(fpr, tpr, color=FOCAL if i == 0 else COMP,
               label=f"{DS_LABEL.get(ds, ds)} ({auc:.2f})")
    c.plot([0, 1], [0, 1], ls=":", color=GREY, lw=0.8)
    c.set_xlim(0, 1); c.set_ylim(0, 1.02)
    c.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); c.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    c.set_xlabel("False positive rate")
    c.set_ylabel("True positive rate")
    c.set_title("Docking as a classifier")
    c.legend(frameon=False, loc="lower right")

    # d - enrichment vs ceiling
    ds_keys = [k for k in P["by_dataset"]]
    x = np.arange(len(ds_keys))
    w = 0.28
    for j, (key, lab, col) in enumerate([("ef5", "EF 5 %", FOCAL),
                                         ("ef10", "EF 10 %", LIGHT)]):
        d.bar(x + (j - 0.5) * w, [P["by_dataset"][k][key] for k in ds_keys],
              width=w, color=col, label=lab)
    for i, k in enumerate(ds_keys):
        ceil = P["by_dataset"][k]["ef_max"]
        d.plot([i - 0.55, i + 0.55], [ceil, ceil], color=COMP, lw=1.1, ls="--")
    d.axhline(1.0, color=GREY, ls=":", lw=0.8)
    d.text(len(ds_keys) - 0.5, 1.04, "no enrichment", fontsize=5.4, color=GREY, ha="right")
    d.set_xticks(x)
    d.set_xticklabels([DS_LABEL.get(k, k).replace(" ", "\n") for k in ds_keys], fontsize=6)
    d.set_ylabel("enrichment factor")
    d.set_title("Enrichment at the top of the list")
    d.set_ylim(0, max(d.get_ylim()[1], max(P["by_dataset"][k]["ef_max"]
                                            for k in ds_keys) * 1.22))
    d.legend(frameon=False, loc="upper left", title="dashed = attainable ceiling",
             title_fontsize=5.4)

    # e - head-to-head on held-out molecules
    H = P.get("held_out_head_to_head", {})
    usable = {k: v for k, v in H.items() if "docking" in v}
    if usable:
        keys = list(usable)
        x = np.arange(len(keys))
        for j, (meth, lab, col) in enumerate([("docking", "docking score", FOCAL),
                                              ("qsar", "QSAR (random forest)", COMP)]):
            e.bar(x + (j - 0.5) * w, [usable[k][meth]["roc_auc"] for k in keys],
                  width=w, color=col, label=lab)
        e.axhline(0.5, color=GREY, ls=":", lw=0.8)
        e.set_xticks(x)
        e.set_xticklabels([f"{DS_LABEL.get(k, k)}\n(n={usable[k]['docking']['n']})"
                           for k in keys], fontsize=6)
        e.set_ylim(0, 1.05)
        e.set_ylabel("ROC-AUC")
        e.legend(frameon=False, loc="upper right")
    else:
        e.text(0.5, 0.5, "too few held-out molecules\nfor a fair comparison",
               transform=e.transAxes, ha="center", va="center", fontsize=6.4, color=GREY)
        e.set_xticks([]); e.set_yticks([])
    e.set_title("Same molecules, neither fitted")

    # f - size bias
    f.scatter(df.MW, df.vina_kcal, s=7, alpha=0.45, color=FOCAL, edgecolors="none")
    z = np.polyfit(df.MW, df.vina_kcal, 1)
    xs = np.linspace(df.MW.min(), df.MW.max(), 50)
    f.plot(xs, np.polyval(z, xs), color=COMP, lw=1.2)
    r = df[["MW", "vina_kcal"]].corr(method="spearman").iloc[0, 1]
    f.text(0.97, 0.06, f"Spearman ρ = {r:.2f}", transform=f.transAxes,
           ha="right", fontsize=6.2)
    lo, hi = df.MW.min() - 15, df.MW.max() + 15
    f.set_xlim(lo, hi)
    f.set_xticks([t for t in range(100, 1001, 100) if lo <= t <= hi])
    f.set_xlabel("molecular weight (Da)")
    f.set_ylabel("Vina score (kcal/mol)")
    f.set_title("Score rewards size")

    for ax, L in zip([a, b, c, d, e, f], "abcdef"):
        panel_letter(ax, L)
    os.makedirs(f"{D}/figures", exist_ok=True)
    out = f"{D}/figures/fig3_docking.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {out}")
    return fig


if __name__ == "__main__":
    main()
