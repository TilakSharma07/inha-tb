"""Step 4 - train XGBoost / random forest / logistic-regression baseline on both split types.

The scaffold-vs-random comparison is the deliverable: it quantifies how much a random split
overstates generalisation on this data. Every headline number in the README is the
scaffold-split number.

Also writes the applicability-domain table: test compounds binned by maximum Tanimoto
similarity to any training compound, scored within each bin. Accuracy is a function of
chemical distance, not a fixed property of the model.
"""
import json, os
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                             balanced_accuracy_score, roc_curve, precision_recall_curve)
from xgboost import XGBClassifier

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

MODELS = {
    "XGBoost": lambda: XGBClassifier(
        n_estimators=400, max_depth=6, learning_rate=0.05, subsample=0.8,
        colsample_bytree=0.6, reg_lambda=1.5, min_child_weight=3,
        eval_metric="logloss", n_jobs=8, random_state=42),
    "RandomForest": lambda: RandomForestClassifier(
        n_estimators=500, max_depth=None, min_samples_leaf=2, max_features="sqrt",
        n_jobs=8, random_state=42, class_weight="balanced"),
    "LogReg(baseline)": lambda: LogisticRegression(
        C=0.1, max_iter=3000, class_weight="balanced"),
}
AD_BINS = [(0.00, 0.25), (0.25, 0.35), (0.35, 0.50), (0.50, 1.01)]


def enrichment(y, p, frac):
    """Enrichment factor at the top `frac` of the ranked list. Ceiling is 1/active_rate."""
    k = max(1, int(len(y) * frac))
    top = np.argsort(-p)[:k]
    return float(y[top].mean() / y.mean()) if y.mean() > 0 else np.nan


def tanimoto_max(Q, Ref):
    """For each row of Q, the maximum Tanimoto similarity to any row of Ref (binary fps)."""
    Q, Ref = Q.astype(np.float32), Ref.astype(np.float32)
    inter = Q @ Ref.T
    union = Q.sum(1)[:, None] + Ref.sum(1)[None, :] - inter
    return np.nan_to_num(inter / np.maximum(union, 1e-9)).max(1)


def main():
    F = pd.read_csv(f"{D}/features.csv")
    X_fp = np.load(f"{D}/fp_matrix.npy")
    S = json.load(open(f"{D}/splits.json"))
    DCOLS = S["descriptor_cols"]

    rows, curves, ad_rows, ad_raw, preds = [], {}, [], {}, []
    for ds, sp in S["splits"].items():
        gi = np.array(sp["global_idx"])
        y = np.array(sp["y"])
        Xd = F.iloc[gi][DCOLS].to_numpy(float)
        X = np.hstack([X_fp[gi], Xd])
        sm = F.iloc[gi]["smiles"].to_numpy()

        for split in ["scaf", "rand"]:
            tr, te = np.array(sp[f"{split}_train"]), np.array(sp[f"{split}_test"])
            sc = StandardScaler().fit(X[tr])
            for name, mk in MODELS.items():
                m, use_sc = mk(), "LogReg" in name
                m.fit(sc.transform(X[tr]) if use_sc else X[tr], y[tr])
                p = m.predict_proba(sc.transform(X[te]) if use_sc else X[te])[:, 1]
                yh = (p >= 0.5).astype(int)
                base = y[te].mean()
                rows.append(dict(
                    dataset=ds, split=split, model=name, n_train=len(tr), n_test=len(te),
                    test_actives=int(y[te].sum()),
                    ROC_AUC=roc_auc_score(y[te], p),
                    PR_AUC=average_precision_score(y[te], p),
                    MCC=matthews_corrcoef(y[te], yh),
                    BalAcc=balanced_accuracy_score(y[te], yh),
                    EF5=enrichment(y[te], p, .05), EF10=enrichment(y[te], p, .10),
                    EF_max=1 / base if base > 0 else np.nan))
                preds.extend(dict(dataset=ds, split=split, model=name,
                                  smiles=sm[i], y=int(y[te][j]), p_active=float(p[j]))
                             for j, i in enumerate(te))
                if name == "XGBoost":
                    fpr, tpr, _ = roc_curve(y[te], p)
                    pr, rc, _ = precision_recall_curve(y[te], p)
                    curves[f"{ds}|{split}"] = dict(fpr=fpr.tolist(), tpr=tpr.tolist(),
                                                   prec=pr.tolist(), rec=rc.tolist(),
                                                   base=float(base))

        # applicability domain: random forest on the scaffold split
        tr, te = np.array(sp["scaf_train"]), np.array(sp["scaf_test"])
        m = MODELS["RandomForest"]()
        m.fit(X[tr], y[tr])
        p = m.predict_proba(X[te])[:, 1]
        sim = tanimoto_max(X_fp[gi][te], X_fp[gi][tr])
        ad_raw[ds] = dict(sim=sim.tolist(), p=p.tolist(), y=y[te].tolist())
        for lo, hi in AD_BINS:
            k = (sim >= lo) & (sim < hi)
            if k.sum() < 5 or len(set(y[te][k])) < 2:
                continue
            ad_rows.append(dict(dataset=ds, tanimoto_bin=f"{lo:.2f}-{hi:.2f}",
                                n=int(k.sum()), actives=int(y[te][k].sum()),
                                ROC_AUC=roc_auc_score(y[te][k], p[k]),
                                mean_abs_err=float(np.abs(p[k] - y[te][k]).mean())))

    R = pd.DataFrame(rows)
    R.to_csv(f"{D}/metrics.csv", index=False)
    json.dump(curves, open(f"{D}/curves.json", "w"))
    AD = pd.DataFrame(ad_rows)
    AD.to_csv(f"{D}/applicability_domain.csv", index=False)
    json.dump(ad_raw, open(f"{D}/ad_raw.json", "w"))
    P = pd.DataFrame(preds)
    P.to_csv(f"{D}/qsar_predictions.csv", index=False)
    print(f"wrote qsar_predictions.csv ({len(P)} test-set predictions)")

    print(R.round(3).to_string(index=False))
    print("\n--- applicability domain (random forest, scaffold split) ---")
    print(AD.round(3).to_string(index=False))
    piv = R.pivot_table(index=["dataset", "model"], columns="split", values="ROC_AUC")
    print("\n--- random-split inflation (rand - scaf) ---")
    print((piv["rand"] - piv["scaf"]).round(3).to_string())


if __name__ == "__main__":
    main()
