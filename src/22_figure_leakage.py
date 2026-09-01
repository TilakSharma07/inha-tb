"""Step 22 - figure 5: is the generalisation gap leakage, or are singletons just hard?

Draws from data/leakage_decomposition.csv (written by src/21), never from a live re-run,
matching the convention used by src/10 and src/19: figures are drawn from committed tables
so that regenerating a figure cannot silently change what it shows.

Panel a: the structural fact that motivates the question - the scaffold-split test set is
entirely singleton scaffolds, because groups are assigned largest-first to train.
Panel b: ROC-AUC by subset within the RANDOM split, where the training set is held fixed,
with bootstrap CIs. Consistent in direction, and underpowered - the figure shows both.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = Path(__file__).resolve().parent.parent
SRC = D / "data" / "leakage_decomposition.csv"
OUT = D / "figures" / "fig5_leakage.png"
FOCAL, COMP, LIGHT, GREY = "#c0562c", "#2a6b9c", "#9ecae1", "#8a8a8a"
DS_LABEL = {"InhA_enzyme": "InhA enzyme", "Mtb_whole_cell": "$M.\\,tb$ whole-cell"}
GROUPS = [("random_test_has_analogues", "random test: analogues in training", COMP),
          ("random_test_singleton", "random test: singleton scaffold", LIGHT),
          ("scaffold_test", "scaffold test (all singletons)", FOCAL)]


def style():
    plt.rcParams.update({
        "font.size": 7, "axes.titlesize": 7.6, "axes.labelsize": 7,
        "xtick.labelsize": 6.4, "ytick.labelsize": 6.4, "legend.fontsize": 6.2,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    })


def main():
    if not SRC.exists():
        sys.exit(f"missing {SRC.relative_to(D)} - run src/21_leakage_decomposition.py first")
    L = pd.read_csv(SRC)
    style()
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.1))

    # --- panel a: scaffold-group size by partition ---------------------------
    ax = axes[0]
    F = D / "data" / "features.csv"
    S = D / "data" / "splits.json"
    if F.exists() and S.exists():
        feat = pd.read_csv(F)
        sp = json.load(open(S))["splits"]
        labels, data, cols = [], [], []
        for ds in sorted(sp):
            f = feat[feat.dataset == ds].reset_index(drop=True)
            size = f.groupby("scaffold").size()
            for part, nm, c in (("scaf_train", "train", GREY),
                                ("scaf_test", "test", FOCAL)):
                idx = [i for i in sp[ds][part] if i < len(f)]
                labels.append(f"{DS_LABEL.get(ds, ds)}\n{nm}")
                data.append(f.loc[idx, "scaffold"].map(size).dropna().to_numpy())
                cols.append(c)
        pos = np.arange(len(data))
        bp = ax.boxplot(data, positions=pos, widths=0.62, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=1.0),
                        whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7),
                        boxprops=dict(lw=0.7))
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c); patch.set_alpha(0.85)
        for x, d, lab in zip(pos, data, labels):
            ax.plot([x], [d.max()], marker="_", ms=7, color="black", lw=0.8)
            if d.max() == 1:            # whole partition is singletons -> box degenerates
                ax.annotate("all 1", xy=(x, 1), xytext=(x, 1.5), fontsize=5.8,
                            color=FOCAL, ha="center", va="bottom",
                            arrowprops=dict(arrowstyle="-", lw=0.6, color=FOCAL,
                                            shrinkA=0.5, shrinkB=1.5))
        ax.set_xticks(pos)
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_ylabel("molecules sharing the scaffold")
        ax.set_yscale("log")
        ax.set_yticks([1, 3, 10, 30])
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_title("Every scaffold-split test molecule is a singleton",
                     loc="left", pad=6)
        ax.axhline(1, color=GREY, ls=":", lw=0.8, zorder=0)
        ax.text(0.02, 0.97, "$-$  largest group in the partition\n"
                "1  =  no scaffold-mate anywhere", fontsize=5.8,
                color="#333", transform=ax.transAxes, va="top", ha="left",
                linespacing=1.5)
    ax.text(-0.10, 1.05, "a", transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")

    # --- panel b: AUC by subset, inside the random split --------------------
    ax = axes[1]
    mo = "RandomForest"
    sub = L[L.model == mo]
    dss = sorted(sub.dataset.unique())
    w = 0.26
    for j, (key, lab, col) in enumerate(GROUPS):
        xs, ys, los, his = [], [], [], []
        for i, ds in enumerate(dss):
            r = sub[(sub.dataset == ds) & (sub.group == key)]
            if r.empty:
                continue
            xs.append(i + (j - 1) * w)
            ys.append(float(r.roc_auc.iloc[0]))
            los.append(float(r.roc_auc.iloc[0]) - float(r.ci_lo.iloc[0]))
            his.append(float(r.ci_hi.iloc[0]) - float(r.roc_auc.iloc[0]))
        ax.bar(xs, ys, width=w * 0.9, color=col, label=lab, zorder=3)
        ax.errorbar(xs, ys, yerr=[los, his], fmt="none", ecolor="#333", elinewidth=0.8,
                    capsize=2.0, capthick=0.8, zorder=4)
        for x, y in zip(xs, ys):
            ax.annotate(f"{y:.2f}", xy=(x, 0.03), fontsize=5.8, ha="center",
                        va="bottom", color="white", zorder=5)
    ax.axhline(0.5, color=GREY, ls=":", lw=0.8, zorder=1)
    ax.text(0.995, 0.5, "random ranking", fontsize=5.8, color=GREY,
            transform=ax.get_yaxis_transform(), ha="right", va="center",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
    ticklabels = []
    for ds in dss:
        ns = [sub[(sub.dataset == ds) & (sub.group == k)].n for k, _, _ in GROUPS]
        ns = [int(x.iloc[0]) for x in ns if len(x)]
        ticklabels.append(f"{DS_LABEL.get(ds, ds)}\nn = " + " / ".join(str(x) for x in ns))
    ax.set_xticks(np.arange(len(dss)))
    ax.set_xticklabels(ticklabels)
    ax.set_ylim(0, 1.08)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel(f"ROC-AUC ({mo}), 95% bootstrap CI")
    ax.set_title("Same model, same training set", loc="left", pad=6)
    ax.legend(frameon=False, loc="upper center", ncol=1, handlelength=0.9,
              handletextpad=0.45, borderpad=0.1, labelspacing=0.3, fontsize=6.0,
              bbox_to_anchor=(0.5, -0.16))
    ax.text(-0.10, 1.05, "b", transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")

    OUT.parent.mkdir(exist_ok=True)
    fig.tight_layout(pad=0.6)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    g = L.dropna(subset=["delta_leakage"])
    print(f"wrote {OUT.relative_to(D)} "
          f"(leakage median {g.delta_leakage.median():+.3f}, "
          f"{(g.delta_leakage > 0).sum()}/{len(g)} positive)")
    return fig


if __name__ == "__main__":
    main()
