"""Step 19 - figure 4, drawn from data/reproducibility_drift.csv (written by src/18).

The repo's convention is that figures come from saved tables, not from a live re-run, so
this script never re-trains. Run src/18_check_reproducibility.py first.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = Path(__file__).resolve().parent.parent
SRC = D / "data" / "reproducibility_drift.csv"
ATTRIB = D / "data" / "version_attribution.csv"
OUT = D / "figures" / "fig4_reproducibility.png"

MODELS = ["LogReg(baseline)", "RandomForest", "XGBoost"]
LBL = {"LogReg(baseline)": "logistic reg.", "RandomForest": "random forest",
       "XGBoost": "XGBoost"}
DS = {"InhA_enzyme": "InhA enzyme", "Mtb_whole_cell": "whole-cell"}
METRICS = [("ROC_AUC", "ROC-AUC"), ("PR_AUC", "PR-AUC"),
           ("BalAcc", "bal. acc."), ("MCC", "MCC")]
COL_C, COL_R = "#1f4e79", "#c26a1f"
COLS = {"LogReg(baseline)": "#9ecae1", "RandomForest": "#5b9bd5", "XGBoost": COL_R}
GATE = 0.02          # ranking-metric gate in src/18


def main():
    if not SRC.exists():
        sys.exit(f"{SRC.relative_to(D)} missing - run src/18_check_reproducibility.py first")
    m = pd.read_csv(SRC)

    plt.rcParams.update({
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.05))

    # panel a - the generalisation gap, committed vs re-run
    ax = axes[0]
    def gap(suffix):
        p = m.pivot_table(index=["dataset", "model"], columns="split",
                          values=f"ROC_AUC_{suffix}")
        return p["rand"] - p["scaf"]
    gc, gr = gap("ref"), gap("now")
    idx = [(d, mo) for d in DS for mo in MODELS]
    y = np.arange(len(idx))[::-1]
    for yy, k in zip(y, idx):
        if abs(gc[k] - gr[k]) > 1e-9:
            ax.plot([gc[k], gr[k]], [yy, yy], color="#9a9a9a", lw=1.1, zorder=1)
    ax.scatter([gc[k] for k in idx], y, s=64, facecolors="none", edgecolors=COL_C,
               linewidths=1.3, zorder=3, label="as committed")
    ax.scatter([gr[k] for k in idx], y, s=22, color=COL_R, marker="D", zorder=4,
               label="re-run here")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{LBL[mo]}\n{DS[d]}" for d, mo in idx])
    ax.set_xlabel("ROC-AUC inflation (random \u2212 scaffold split)\n"
                  "higher = the random split overstated more", labelpad=11)
    ax.set_title("The generalisation gap survives the re-run", loc="left", pad=6)
    ax.axvline(0, color="#555", lw=0.8, ls=(0, (4, 3)), zorder=0)
    ax.set_xlim(-0.012, 0.142); ax.set_xticks([0.00, 0.04, 0.08, 0.12]); ax.margins(y=0.10)
    ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.02, -0.02),
              handletextpad=0.3, borderpad=0.1, labelspacing=0.3)

    # panel b - drift attributed to each factor. The point: the endpoint-to-endpoint
    # comparison a reviewer would make is NOT the worst case, because the two factors
    # partly cancel. A single re-run therefore does not bound the drift.
    ax = axes[1]
    if not ATTRIB.exists():
        ax.text(0.5, 0.5, "run src/20_version_attribution.py", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="#8c2f2f")
        worst = float("nan")
    else:
        A = pd.read_csv(ATTRIB)
        envs = list(dict.fromkeys(A.env))
        key = ["dataset", "split", "model"]

        def drift(ea, eb):
            x = A[(A.env == ea) & (A.model == "XGBoost")]
            y = A[(A.env == eb) & (A.model == "XGBoost")]
            mm = x.merge(y, on=key, suffixes=("_a", "_b"))
            return [(mm[f"{c}_a"] - mm[f"{c}_b"]).abs().max() for c, _ in METRICS]

        # ordered so the endpoint pair (what a re-run measures) is last
        pairs = [(envs[0], envs[1], "python only"),
                 (envs[1], envs[2], "xgboost+numpy only"),
                 (envs[0], envs[2], "both (a re-run)")]
        pc = ["#9ecae1", "#5b9bd5", COL_R]
        w = 0.26
        for j, ((ea, eb, lab), col) in enumerate(zip(pairs, pc)):
            d = drift(ea, eb)
            ax.bar(np.arange(len(METRICS)) + (j - 1) * w, d, width=w * 0.9,
                   color=col, label=lab, zorder=3)
        worst = max(max(drift(a, b)[2:4]) for a, b, _ in pairs)
        endpoint = drift(envs[0], envs[2])[2]
        mid = drift(envs[0], envs[1])[2]
        ax.annotate("", xy=(2 + w, endpoint), xytext=(2 - w, mid),
                    arrowprops=dict(arrowstyle="->", lw=0.8, color="#444",
                                    shrinkA=1.5, shrinkB=1.5,
                                    connectionstyle="arc3,rad=-0.25"), zorder=6)
        ax.text(2.0, mid * 1.06, "a re-run understates it", fontsize=6, color="#444",
                ha="center", va="bottom")

    ax.set_xticks(np.arange(len(METRICS)))
    ax.set_xticklabels([l for _, l in METRICS])
    ax.set_ylabel("largest XGBoost drift", labelpad=4)
    ax.set_title("Version effects partly cancel", loc="left", pad=6)
    ax.set_ylim(0, 0.058); ax.set_yticks([0, 0.02, 0.04])
    # The 0.02 gate applies ONLY to the ranking metrics; MCC and balanced accuracy are
    # reported against a looser 0.10. Spanning the whole panel would imply the two
    # right-hand groups are in breach, which they are not.
    ax.plot([-0.5, 1.5], [GATE, GATE], color="#8c2f2f", lw=0.9, ls=(0, (3, 2)), zorder=2,
            clip_on=False)
    ax.text(-0.45, GATE * 1.04, "gate (ranking metrics only)", fontsize=6,
            color="#8c2f2f", va="bottom")
    ax.text(2.5, 0.0555, "reported, gated at 0.10", fontsize=6, color="#777",
            ha="center", va="top")
    ax.plot([1.5, 1.5], [0, 0.052], color="#ccc", lw=0.7, zorder=1)
    ax.legend(frameon=False, loc="upper left", handlelength=0.9, handletextpad=0.45,
              borderpad=0.15, labelspacing=0.3, title="varied:", title_fontsize=6)

    for ax_, letter in zip(axes, "ab"):
        ax_.text(-0.02, 1.06, letter, transform=ax_.transAxes, fontsize=10,
                 fontweight="bold", va="bottom", ha="right")

    fig.tight_layout(w_pad=2.2)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(D)} (worst threshold-metric drift {worst:.4f})")


if __name__ == "__main__":
    main()
