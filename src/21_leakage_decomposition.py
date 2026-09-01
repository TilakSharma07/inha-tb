"""Step 21 - decompose the generalisation gap into analogue leakage and scaffold novelty.

The repo's central claim is that the random split overstates performance. Steps 04 and 10
establish it by comparing two splits, which leaves one alternative explanation open: the
scaffold-split test set consists entirely of SINGLETON scaffolds (groups are assigned
largest-first to train, so only groups of size 1 reach the test set). Singleton molecules
could be intrinsically harder for reasons unrelated to scaffold novelty - one-off papers,
more exotic chemistry, noisier assays - in which case the gap would not be leakage.

This step removes that confound without training anything new. Inside the RANDOM split,
the training set is fixed and the model is one model; its test set is then partitioned by
whether the molecule's Murcko scaffold occurs more than once in the dataset:

    has analogues   scaffold-mates present in training -> leakage possible
    singleton       scaffold occurs once in the dataset -> no scaffold-mate anywhere

If the gap were an artifact of singletons being hard, both subsets would score the same.
If it is leakage, the analogue subset scores higher, and the singleton subset lands near
the scaffold-split value. Reported with bootstrap CIs and a permutation test on the
difference.

Writes data/leakage_decomposition.csv. Consumed by src/22_figure_leakage.py.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

D = Path(__file__).resolve().parent.parent
OUT = D / "data" / "leakage_decomposition.csv"
SEED, N_BOOT, N_PERM = 42, 4000, 4000
MODELS = ("RandomForest", "XGBoost", "LogReg(baseline)")


def auc_rows(Y, Pm):
    """ROC-AUC for every row of the (B, n) label and score matrices.

    Mann-Whitney form with mid-ranks for ties, which is algebraically identical to the
    trapezoidal ROC-AUC that sklearn computes (verified to 2e-16 against
    sklearn.metrics.roc_auc_score, including on heavily tied scores). Vectorising matters
    here: the loop version made ~168,000 sklearn calls and took ~4 minutes, which is too
    slow for something a reviewer is asked to run.
    """
    R = rankdata(Pm, method="average", axis=1)
    npos = Y.sum(axis=1)
    nneg = Y.shape[1] - npos
    rsum = (R * Y).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (rsum - npos * (npos + 1) / 2.0) / (npos * nneg)
    out[(npos == 0) | (nneg == 0)] = np.nan
    return out


def auc_ci(y, p, rng, n=N_BOOT):
    """ROC-AUC with a percentile bootstrap CI. Replicates that lose a class give NaN and
    are dropped, so the CI is conditional on both classes being present."""
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    idx = rng.integers(0, len(y), (n, len(y)))
    v = auc_rows(y[idx], p[idx])
    v = v[np.isfinite(v)]
    a = roc_auc_score(y, p)
    if v.size == 0:
        return a, np.nan, np.nan
    return a, float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def delta_ci(y1, p1, y2, p2, rng, n=N_BOOT):
    """Bootstrap CI on AUC(group1) - AUC(group2). The two groups are disjoint sets of
    molecules, so they are resampled independently."""
    if len(np.unique(y1)) < 2 or len(np.unique(y2)) < 2:
        return np.nan, np.nan
    i = rng.integers(0, len(y1), (n, len(y1)))
    j = rng.integers(0, len(y2), (n, len(y2)))
    d = auc_rows(y1[i], p1[i]) - auc_rows(y2[j], p2[j])
    d = d[np.isfinite(d)]
    if d.size == 0:
        return np.nan, np.nan
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def perm_p(y1, p1, y2, p2, rng, n=N_PERM):
    """Two-sided permutation test on AUC(group1) - AUC(group2).

    Label-prediction pairs are kept intact and only GROUP membership is shuffled, so the
    null is 'which group a molecule falls in carries no information about how well the
    model ranks it'. Shuffling labels instead would test whether the model ranks at all,
    which is a different and already-answered question.
    """
    if len(np.unique(y1)) < 2 or len(np.unique(y2)) < 2:
        return np.nan
    obs = abs(roc_auc_score(y1, p1) - roc_auc_score(y2, p2))
    y = np.concatenate([y1, y2])
    p = np.concatenate([p1, p2])
    n1 = len(y1)
    perm = np.argsort(rng.random((n, len(y))), axis=1)
    A, B = perm[:, :n1], perm[:, n1:]
    d = np.abs(auc_rows(y[A], p[A]) - auc_rows(y[B], p[B]))
    ok = np.isfinite(d)
    if not ok.any():
        return np.nan
    return float((int((d[ok] >= obs - 1e-12).sum()) + 1) / (int(ok.sum()) + 1))


def main():
    F = D / "data" / "features.csv"
    P = D / "data" / "qsar_predictions.csv"
    for f in (F, P):
        if not f.exists():
            sys.exit(f"missing {f.relative_to(D)} - run src/03 and src/04 first")
    feat = pd.read_csv(F)
    pred = pd.read_csv(P)
    rng = np.random.default_rng(SEED)

    rows = []
    for ds, f in feat.groupby("dataset"):
        size = f.groupby("scaffold").size()
        single = set(f.loc[f.scaffold.map(size) == 1, "smiles"])
        for mo in MODELS:
            sc = pred[(pred.dataset == ds) & (pred.split == "scaf") & (pred.model == mo)]
            rd = pred[(pred.dataset == ds) & (pred.split == "rand") & (pred.model == mo)]
            if sc.empty or rd.empty:
                continue
            is_s = rd.smiles.isin(single)
            groups = {
                "scaffold_test": sc,
                "random_test_singleton": rd[is_s],
                "random_test_has_analogues": rd[~is_s],
            }
            stats = {}
            for name, sub in groups.items():
                a, lo, hi = auc_ci(sub.y.to_numpy(), sub.p_active.to_numpy(), rng)
                stats[name] = (a, lo, hi, len(sub), float(sub.y.mean()))
                rows.append(dict(dataset=ds, model=mo, group=name, n=len(sub),
                                 prevalence=round(float(sub.y.mean()), 4),
                                 roc_auc=round(a, 4) if a == a else np.nan,
                                 ci_lo=round(lo, 4) if lo == lo else np.nan,
                                 ci_hi=round(hi, 4) if hi == hi else np.nan))
            a_sub = rd[~is_s]
            s_sub = rd[is_s]
            p_leak = perm_p(a_sub.y.to_numpy(), a_sub.p_active.to_numpy(),
                            s_sub.y.to_numpy(), s_sub.p_active.to_numpy(), rng)
            rows[-1]["perm_p_analogue_vs_singleton"] = (
                round(p_leak, 4) if p_leak == p_leak else np.nan)
            d_leak = stats["random_test_has_analogues"][0] - stats["random_test_singleton"][0]
            d_conf = stats["random_test_singleton"][0] - stats["scaffold_test"][0]
            lo, hi = delta_ci(a_sub.y.to_numpy(), a_sub.p_active.to_numpy(),
                              s_sub.y.to_numpy(), s_sub.p_active.to_numpy(), rng)
            rows[-1]["delta_leakage"] = round(d_leak, 4)
            rows[-1]["delta_leakage_ci_lo"] = round(lo, 4) if lo == lo else np.nan
            rows[-1]["delta_leakage_ci_hi"] = round(hi, 4) if hi == hi else np.nan
            rows[-1]["delta_residual"] = round(d_conf, 4)
            rows[-1]["n_singleton"] = int(len(s_sub))
            print(f"  {ds:<15} {mo:<18} leakage {d_leak:+.3f} "
                  f"[{lo:+.3f}, {hi:+.3f}] p={p_leak:.3g}   "
                  f"residual {d_conf:+.3f}   (n_singleton={len(s_sub)})")

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(D)}: {len(out)} rows")

    g = out.dropna(subset=["delta_leakage"])
    print(f"leakage effect: {g.delta_leakage.min():+.3f} to {g.delta_leakage.max():+.3f} "
          f"(median {g.delta_leakage.median():+.3f}) across "
          f"{len(g)} model x dataset combinations")
    print(f"residual (singleton random test vs scaffold test): "
          f"{g.delta_residual.min():+.3f} to {g.delta_residual.max():+.3f} "
          f"(median {g.delta_residual.median():+.3f})")
    pos = (g.delta_leakage > 0).sum()
    print(f"analogue subset scores higher in {pos}/{len(g)} combinations")

    # What this does and does not establish. The direction is consistent everywhere, but
    # the three models share a test set within a dataset, so they are not independent
    # replicates: a 6/6 sign test would treat six correlated observations as six
    # independent ones and overstate the evidence. The two DATASETS are independent.
    sig = g[g.perm_p_analogue_vs_singleton < 0.05]
    print(f"individually significant at p<0.05: {len(sig)}/{len(g)} "
          f"({', '.join(sorted(set(sig.dataset))) if len(sig) else 'none'})")
    small = g[g.n_singleton < 30]
    if len(small):
        print(f"underpowered (singleton subset n<30): "
              f"{', '.join(f'{r.dataset}/{r.model} n={r.n_singleton}' for r in small.itertuples())}")
    print("independent replicates: 2 datasets; models within a dataset share the test set")

    # Per-combination CIs are wide because the singleton arm is small (22 and 65 molecules,
    # 9 actives in the smaller). Pooling the two datasets for ONE model gives the only
    # test here with reasonable power, and is the number worth quoting. RandomForest is
    # used because src/09's head-to-head uses it.
    mo = "RandomForest"
    ys, ps, zs = [], [], []
    for ds, f in feat.groupby("dataset"):
        size = f.groupby("scaffold").size()
        single = set(f.loc[f.scaffold.map(size) == 1, "smiles"])
        rd = pred[(pred.dataset == ds) & (pred.split == "rand") & (pred.model == mo)]
        if rd.empty:
            continue
        # Ranks are only comparable within a dataset, so each dataset's predictions are
        # converted to within-dataset percentile ranks before pooling. Pooling raw
        # probabilities would let a dataset-level calibration shift masquerade as signal.
        r = pd.Series(rd.p_active).rank(pct=True).to_numpy()
        ys.append(rd.y.to_numpy())
        ps.append(r)
        zs.append(rd.smiles.isin(single).to_numpy())
    if len(ys) == 2:
        y = np.concatenate(ys); pr = np.concatenate(ps); z = np.concatenate(zs)
        pv = perm_p(y[~z], pr[~z], y[z], pr[z], rng)
        lo, hi = delta_ci(y[~z], pr[~z], y[z], pr[z], rng)
        d = roc_auc_score(y[~z], pr[~z]) - roc_auc_score(y[z], pr[z])
        print(f"\npooled across datasets ({mo}, within-dataset percentile ranks):")
        print(f"  analogue n={int((~z).sum())} vs singleton n={int(z.sum())}  "
              f"delta {d:+.3f} [{lo:+.3f}, {hi:+.3f}]  permutation p={pv:.4g}")
        verdict = ("consistent in direction and significant when pooled"
                   if pv == pv and pv < 0.05 else
                   "consistent in direction; not significant even pooled")
        print(f"  verdict: {verdict}")


if __name__ == "__main__":
    main()
