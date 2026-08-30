"""Step 9 - join docking scores to measured potency; compare docking against the QSAR model.

This is where most screening write-ups overclaim. Vina's empirical scoring function was
fitted to reproduce binding modes, not affinities, so the honest questions are narrow:

  1. Does the score separate measured actives from inactives at all (ROC-AUC)?
  2. Does it enrich actives at the top of the list (EF at 5 % / 10 %)?
  3. How does that compare to the ligand-based QSAR model on the same molecules?

Question 3 is the one that makes the answer interpretable, and it has a trap: the QSAR model
saw most of these molecules during training, so a global comparison flatters QSAR. The
comparison here is therefore restricted to molecules in the QSAR *scaffold-split test set* —
molecules neither method has fitted. That subset is small, so its numbers are reported with
their n and treated as indicative rather than conclusive.

A weak-but-positive docking result is the expected outcome and is reported as such.
"""
import json, os
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CUTOFF = 6.0


def enrichment(y, score_desc, frac):
    """EF at top `frac`. score_desc: higher = predicted better."""
    k = max(1, int(len(y) * frac))
    top = np.argsort(-score_desc)[:k]
    return float(y[top].mean() / y.mean()) if y.mean() > 0 else np.nan


def perf(y, score_desc):
    return dict(n=int(len(y)), actives=int(y.sum()),
                roc_auc=float(roc_auc_score(y, score_desc)),
                ef5=enrichment(y, score_desc, 0.05),
                ef10=enrichment(y, score_desc, 0.10),
                ef_max=float(len(y) / y.sum()))


def main():
    src = f"{D}/dock/docking_scores.jsonl"
    rows = ([json.loads(l) for l in open(src)] if os.path.exists(src)
            else json.load(open(f"{D}/dock/docking_scores.json")))
    S = pd.DataFrame(rows)
    S = S[S.vina_kcal.notna()].drop_duplicates("lig_id").copy()
    lib = pd.read_csv(f"{D}/dock/library.csv")
    df = lib.merge(S, on="lig_id", how="inner")
    df["active"] = (df.pchembl >= CUTOFF).astype(int)
    df["vina_desc"] = -df.vina_kcal            # more negative kcal = better -> flip sign
    print(f"{len(df)}/{len(lib)} library molecules scored", flush=True)

    out = {}
    for ds, sub in df.groupby("dataset"):
        y = sub.active.to_numpy()
        if y.sum() < 5 or y.sum() == len(y):
            print(f"{ds}: skipped (actives={y.sum()}/{len(y)})")
            continue
        out[ds] = perf(y, sub.vina_desc.to_numpy())
        out[ds].update(
            spearman=float(sub[["vina_kcal", "pchembl"]].corr(method="spearman").iloc[0, 1]),
            pearson=float(sub[["vina_kcal", "pchembl"]].corr().iloc[0, 1]),
            vina_median=float(sub.vina_kcal.median()),
            vina_best=float(sub.vina_kcal.min()))
        print(f"{ds}: n={out[ds]['n']} act={out[ds]['actives']} "
              f"ROC-AUC={out[ds]['roc_auc']:.3f} rho={out[ds]['spearman']:.3f} "
              f"EF5={out[ds]['ef5']:.2f} (ceiling {out[ds]['ef_max']:.2f})")

    # --- head-to-head on molecules neither method fitted -------------------------------
    head = {}
    pf = f"{D}/data/qsar_predictions.csv"
    if os.path.exists(pf):
        P = pd.read_csv(pf)
        P = P[(P.split == "scaf") & (P.model == "RandomForest")]
        for ds, sub in df.groupby("dataset"):
            pj = P[P.dataset == ds][["smiles", "p_active"]]
            m = sub.merge(pj, on="smiles", how="inner")
            y = m.active.to_numpy()
            if len(m) < 20 or y.sum() < 5 or y.sum() == len(y):
                head[ds] = dict(n=int(len(m)), actives=int(y.sum()),
                                note="too few held-out molecules for a meaningful comparison")
                print(f"{ds}: head-to-head skipped (n={len(m)}, actives={y.sum()})")
                continue
            head[ds] = dict(docking=perf(y, m.vina_desc.to_numpy()),
                            qsar=perf(y, m.p_active.to_numpy()))
            print(f"{ds}: held-out n={len(m)} | docking ROC-AUC "
                  f"{head[ds]['docking']['roc_auc']:.3f} vs QSAR "
                  f"{head[ds]['qsar']['roc_auc']:.3f}")
    else:
        print("qsar_predictions.csv absent - run src/04_train_qsar.py for the comparison")

    df.to_csv(f"{D}/data/docking_joined.csv", index=False)
    json.dump(dict(by_dataset=out, held_out_head_to_head=head, cutoff=CUTOFF,
                   n_scored=int(len(df)), n_library=int(len(lib))),
              open(f"{D}/data/docking_performance.json", "w"), indent=1)


if __name__ == "__main__":
    main()
