"""Step 18 - re-train from the committed data and report drift against the committed metrics.

Why this exists: a reviewer who re-runs this pipeline should be able to tell whether the
numbers they get are the numbers the repo reports. They will not be bit-identical, and the
honest thing is to say so with a bound rather than to assert equality that does not hold.

What actually happens across library versions (measured, not assumed):

  * splits                      invariant - the Bemis-Murcko split is deterministic by
                                construction (sorted by group size, then name), not by seed
  * RandomForest, logistic reg. bit-identical
  * XGBoost ranking metrics     up to ~0.01 ROC-AUC / PR-AUC between xgboost 3.2.0 and
                                3.4.1 - a MINOR version bump, not a major one
  * XGBoost MCC, BalAcc         up to ~0.03; these threshold the probability at 0.5, so a
                                small shift relabels molecules sitting near the boundary
  * EF5                         coarsely quantised, NOT gated - at n_test = 88 the top 5%
                                is 4 molecules, so one rank change moves EF5 by ~0.35

Gating policy: this exits non-zero if the split changes, if a ranking metric drifts beyond
tolerance, or if the repo's central claim fails - that scaffold-split ROC-AUC is below
random-split ROC-AUC in every model x dataset combination. Metric drift within tolerance is
reported, not treated as failure. The run happens in a temporary copy, so the committed
tables are never touched.
"""
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# The pipeline steps run as subprocesses under this same interpreter. Check before the
# third-party imports below, or a bare ModuleNotFoundError traceback is all the user sees.
NEEDED = ("rdkit", "sklearn", "xgboost", "pandas", "numpy", "scipy")
_missing = [m for m in NEEDED if importlib.util.find_spec(m) is None]
if _missing:
    sys.exit(
        f"this interpreter cannot run the pipeline - missing: {', '.join(_missing)}\n"
        f"  interpreter: {sys.executable}\n"
        f"  fix: conda env create -f environment.yml && conda activate inha-tb\n"
        f"  then re-run: python src/18_check_reproducibility.py"
    )

import pandas as pd

D = Path(__file__).resolve().parent.parent
LOCK = D / "environment.lock.json"

# Tolerances, and where each number comes from. Measured across three interpreters
# (see src/20_version_attribution.py), not chosen to make this check pass:
#
#   ROC_AUC, PR_AUC   worst observed 0.010 / 0.014  -> gate 0.02
#   MCC, BalAcc       worst observed 0.047 / 0.023  -> reported, gated only at 0.10
#   EF5               worst observed 0.75           -> reported, never gated
#
# ROC-AUC and PR-AUC are properties of the ranking, so a code defect shows up in them.
# MCC and balanced accuracy additionally depend on where the predicted probability falls
# relative to a fixed 0.5 cut, so a shift of ~1e-3 in one probability relabels a molecule
# and moves them by ~1/n_test. On an 88-molecule held-out set that is arithmetic noise
# rather than a change in model quality: gating them tightly would be gating the arithmetic.
# The 0.10 gate is deliberately loose - it catches a model that has genuinely changed, not
# version noise. It is set from the worst case across the three interpreters tested, with
# headroom; a wider version sweep could justify moving it.
TOL = {"ROC_AUC": 0.02, "PR_AUC": 0.02, "MCC": 0.10, "BalAcc": 0.10}
GATED = ("ROC_AUC", "PR_AUC")            # a breach here fails the check
REPORTED = ("MCC", "BalAcc")             # printed, gated only against gross change


def reference_status(obs):
    """Compare the running interpreter against the versions that produced the committed
    tables. Being off-reference is not an error - it is the condition this check measures."""
    if not LOCK.exists():
        return None, []
    ref = {k: v for k, v in json.load(open(LOCK)).items() if not k.startswith("_")}
    alias = {"sklearn": "scikit-learn"}
    diff = [(k, ref[alias.get(k, k)], v) for k, v in obs.items()
            if alias.get(k, k) in ref and ref[alias.get(k, k)] != v]
    return ref, sorted(diff)


def versions():
    import numpy, scipy, sklearn, xgboost
    return {"python": sys.version.split()[0], "numpy": numpy.__version__,
            "scipy": scipy.__version__, "sklearn": sklearn.__version__,
            "xgboost": xgboost.__version__}


def committed(relpath):
    """Return the committed (HEAD) bytes of a tracked file, or None.

    The reference for this check has to be the metrics that were committed, not the
    ones sitting in the working tree: run_all.sh runs step 04 before step 18, so step
    04 overwrites data/metrics.csv with output from THIS interpreter. Reading the
    working tree there makes the check compare the re-run against itself - it reports
    'identical' for every model and can never fail. Measured on xgboost 3.2.0 -> 3.4.1:
    working-tree reference reported 0.000 drift where the committed reference reported
    ROC_AUC 0.0099 / EF5 0.376.
    """
    try:
        r = subprocess.run(["git", "show", f"HEAD:{relpath}"], cwd=str(D),
                           capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout if r.returncode == 0 and r.stdout else None


def main():
    import io

    # Reference = the committed tables. Fall back to the working tree only if this is
    # not a git checkout, and say so, because that fallback weakens the check.
    ref_src = "HEAD"
    blob = committed("data/metrics.csv")
    if blob is not None:
        ref = pd.read_csv(io.BytesIO(blob))
    else:
        ref_path = D / "data" / "metrics.csv"
        if not ref_path.exists():
            sys.exit("data/metrics.csv missing - run src/04_train_qsar.py first")
        ref = pd.read_csv(ref_path)
        ref_src = "working tree"

    sp_blob = committed("data/splits.json")
    if sp_blob is not None:
        split_ref = json.loads(sp_blob)["splits"]
    else:
        sp_ref = D / "data" / "splits.json"
        split_ref = json.load(open(sp_ref))["splits"] if sp_ref.exists() else None
        if ref_src == "HEAD":
            ref_src = "HEAD (metrics) + working tree (splits)"

    with tempfile.TemporaryDirectory(prefix="inha-repro-") as tmp:
        W = Path(tmp) / "repo"
        # copy only what steps 3 and 4 need, so nothing in the real tree can be written
        (W / "data").mkdir(parents=True)
        shutil.copytree(D / "src", W / "src")
        for f in ("chembl_filtered.csv", "activity_aggregated.csv"):
            src = D / "data" / f
            if src.exists():
                shutil.copy(src, W / "data" / f)

        for step in ("03_features.py", "04_train_qsar.py"):
            r = subprocess.run([sys.executable, str(W / "src" / step)],
                               capture_output=True, text=True, cwd=str(W))
            if r.returncode:
                sys.exit(f"src/{step} failed in the isolated copy:\n{r.stderr[-800:]}")

        now = pd.read_csv(W / "data" / "metrics.csv")
        sp_now = W / "data" / "splits.json"
        split_now = json.load(open(sp_now))["splits"] if sp_now.exists() else None

    key = ["dataset", "split", "model"]
    m = ref.merge(now, on=key, suffixes=("_ref", "_now"))
    # persist the comparison so src/19 can draw it from a table rather than from
    # a re-run - the repo's convention is that figures come from saved tables.
    # Stamp WHICH reference produced it: a table built against the working tree is
    # a comparison of a re-run against itself (all-zero drift), and figure 4 drawn
    # from such a table is a false claim of bit-reproducibility. src/19 refuses it.
    m["reference_source"] = ref_src
    m.to_csv(D / "data" / "reproducibility_drift.csv", index=False)
    if len(m) != len(ref):
        sys.exit(f"row mismatch: {len(m)} matched of {len(ref)} - the model grid changed")

    obs = versions()
    print(f"reproducibility: {json.dumps(obs)}")
    if ref_src == "HEAD":
        print("  reference: committed tables at HEAD (not the working tree, which step 04 "
              "overwrites)")
    else:
        print(f"  WEAKENED: reference read from the {ref_src} - if step 04 ran in this "
              "interpreter the comparison is partly against itself")
    lock, diff = reference_status(obs)
    if lock is None:
        print("  environment.lock.json missing - cannot say whether this is the reference env")
    elif not diff:
        print("  IN REFERENCE: same versions that produced the committed tables; "
              "any drift below is a real defect, not a version effect")
    else:
        print("  OFF REFERENCE: " + ", ".join(f"{k} {r} -> {o}" for k, r, o in diff))
        print("  drift below is attributable to these, not to the code")

    for c in ("n_train", "n_test", "test_actives"):
        d = (m[f"{c}_ref"] - m[f"{c}_now"]).abs().max()
        if d != 0:
            sys.exit(f"FAIL: {c} differs by up to {d} - the split is not reproducible")

    # Equal sizes are not equal splits: two different partitions can have identical counts.
    # Compare membership index-by-index against the committed data/splits.json.
    if split_ref and split_now:
        for ds in sorted(split_ref):
            for k in ("scaf_train", "scaf_test", "rand_train", "rand_test"):
                a, b = split_ref[ds].get(k), split_now[ds].get(k)
                if a is None or b is None:
                    continue
                if a != b:
                    moved = len(set(map(str, a)) ^ set(map(str, b))) // 2
                    order_only = sorted(map(str, a)) == sorted(map(str, b))
                    sys.exit(f"FAIL: {ds} {k} membership changed "
                             f"({moved} molecules moved, order-only={order_only})")
        print("  splits invariant (sizes AND membership identical, all four partitions)")
    else:
        print("  splits invariant (sizes only - data/splits.json missing, "
              "membership NOT verified)")

    breaches = []
    for mo in sorted(m.model.unique()):
        s = m[m.model == mo]
        parts, worst = [], 0.0
        for c in GATED + REPORTED:
            d = (s[f"{c}_ref"] - s[f"{c}_now"]).abs().max()
            if c in GATED:
                worst = max(worst, d / TOL[c])
            if d > TOL[c]:
                breaches.append(f"{mo} {c} drift {d:.4f} > tol {TOL[c]}")
            if d:
                parts.append(f"{c} {d:.4f}{'' if c in GATED else '*'}")
        print(f"  {mo:18s} {'identical' if not parts else 'drift: ' + ', '.join(parts)}"
              f"{'' if worst <= 1 else '   BREACH'}")

    ef = (m["EF5_ref"] - m["EF5_now"]).abs().max()
    print(f"  EF5 drift {ef:.3f} (not gated - quantised by 1/(0.05*n_test))")

    # The README states "up to X" bounds. A bound that the measurement exceeds is a false
    # claim, so verify the direction here rather than trusting that it was rounded up.
    attrib = D / "data" / "version_attribution.csv"
    readme = D / "README.md"
    if attrib.exists() and readme.exists():
        import re
        A = pd.read_csv(attrib)
        envs = list(dict.fromkeys(A.env))
        if len(envs) >= 2:
            txt = readme.read_text(encoding="utf-8")
            claims = {}
            mm = re.search(r"ROC-AUC, PR-AUC \| up to (\d+\.\d+) / (\d+\.\d+)", txt)
            if mm:
                claims["ROC_AUC"], claims["PR_AUC"] = float(mm.group(1)), float(mm.group(2))
            mm = re.search(r"MCC, balanced accuracy \| up to (\d+\.\d+) / (\d+\.\d+)", txt)
            if mm:
                claims["MCC"], claims["BalAcc"] = float(mm.group(1)), float(mm.group(2))
            mm = re.search(r"\| EF5 \| up to (\d+\.\d+)", txt)
            if mm:
                claims["EF5"] = float(mm.group(1))
            bad = []
            for c, bound in claims.items():
                w = 0.0
                for i, ea in enumerate(envs):
                    for eb in envs[i + 1:]:
                        x = A[(A.env == ea) & (A.model == "XGBoost")]
                        y = A[(A.env == eb) & (A.model == "XGBoost")]
                        q = x.merge(y, on=key, suffixes=("_a", "_b"))
                        w = max(w, (q[f"{c}_a"] - q[f"{c}_b"]).abs().max())
                if w > bound + 1e-9:
                    bad.append(f"{c}: README says up to {bound}, measured {w:.4f}")
            if bad:
                sys.exit("FAIL: README states bounds the measurement exceeds:\n  "
                         + "\n  ".join(bad))
            print(f"  README bounds hold for all {len(claims)} stated metrics")

    p = now.pivot_table(index=["dataset", "model"], columns="split", values="ROC_AUC")
    gap = p["rand"] - p["scaf"]
    print(f"  random - scaffold gap: {gap.min():.3f} to {gap.max():.3f}, "
          f"median {gap.median():.3f}")
    if (gap <= 0).any():
        sys.exit(f"FAIL: scaffold split not below random for {list(gap[gap <= 0].index)} - "
                 "the generalisation-gap claim does not hold under these versions")
    print(f"  central claim holds in all {len(gap)} model x dataset combinations")

    if breaches:
        sys.exit("FAIL:\n  " + "\n  ".join(breaches))


if __name__ == "__main__":
    main()
