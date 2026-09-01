"""Step 20 - attribute metric drift to specific library versions.

NOT part of run_all.sh: this needs several interpreters with different library versions,
which a reviewer will not have. The result is committed as data/version_attribution.csv so
that src/19 can draw figure 4 without re-running it, and so the README's stated bounds are
traceable to a measurement rather than to a single observation.

Usage:
    python src/20_version_attribution.py <label>=<python> [<label>=<python> ...]

Each interpreter re-runs steps 03 and 04 in its own temporary copy of the repo, so the
committed tables are never touched. What varies is the interpreter; the input data does not.

Why this exists: the first version of the README blamed XGBoost drift on a major-version
bump. Holding xgboost fixed at 3.2.0 and moving only python 3.11 -> 3.13 turned out to move
MCC further (0.047) than the full version jump did (0.028), so the single-observation bound
in the README understated the worst case. Intermediate version combinations are not
bracketed by the endpoints.
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

D = Path(__file__).resolve().parent.parent
OUT = D / "data" / "version_attribution.csv"
STEPS = ("03_features.py", "04_train_qsar.py")
PROBE = ("import json,sys,platform;"
         "import numpy,scipy,sklearn,xgboost,rdkit;"
         "print(json.dumps({'python':platform.python_version(),'numpy':numpy.__version__,"
         "'scipy':scipy.__version__,'sklearn':sklearn.__version__,"
         "'xgboost':xgboost.__version__,'rdkit':rdkit.__version__}))")


def probe(py):
    r = subprocess.run([py, "-c", PROBE], capture_output=True, text=True)
    if r.returncode:
        return None
    return json.loads(r.stdout)


def run_one(label, py, vers):
    import pandas as pd
    with tempfile.TemporaryDirectory(prefix=f"inha-attrib-{label}-") as tmp:
        W = Path(tmp) / "repo"
        (W / "data").mkdir(parents=True)
        shutil.copytree(D / "src", W / "src")
        for f in ("chembl_filtered.csv", "activity_aggregated.csv"):
            if (D / "data" / f).exists():
                shutil.copy(D / "data" / f, W / "data" / f)
        for step in STEPS:
            r = subprocess.run([py, str(W / "src" / step)],
                               capture_output=True, text=True, cwd=str(W))
            if r.returncode:
                print(f"  {label}: {step} failed - {r.stderr.strip().splitlines()[-1][:90]}")
                return None
        m = pd.read_csv(W / "data" / "metrics.csv")
    m.insert(0, "env", label)
    for k, v in vers.items():
        m[f"v_{k}"] = v
    return m


def main(argv):
    import pandas as pd
    if len(argv) < 2:
        sys.exit(__doc__.split("Usage:")[1].split("Each")[0].strip())
    specs = []
    for a in argv:
        if "=" not in a:
            sys.exit(f"expected <label>=<python>, got {a!r}")
        label, py = a.split("=", 1)
        vers = probe(py)
        if vers is None:
            print(f"  skipping {label}: {py} cannot import the pipeline's dependencies")
            continue
        specs.append((label, py, vers))
    if len(specs) < 2:
        sys.exit(f"need at least 2 usable interpreters, got {len(specs)}")

    frames = []
    for label, py, vers in specs:
        print(f"  {label:<10} python {vers['python']:<8} xgboost {vers['xgboost']:<7} "
              f"numpy {vers['numpy']:<7} rdkit {vers['rdkit']}")
        got = run_one(label, py, vers)
        if got is not None:
            frames.append(got)
    if len(frames) < 2:
        sys.exit("fewer than 2 interpreters produced metrics")

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT.relative_to(D)}: {out.env.nunique()} environments x "
          f"{len(out) // out.env.nunique()} rows")

    # the number the README must quote: worst drift between ANY pair of environments
    key = ["dataset", "split", "model"]
    envs = list(out.env.unique())
    print("\nworst drift between any pair of environments (XGBoost):")
    for c in ("ROC_AUC", "PR_AUC", "MCC", "BalAcc", "EF5"):
        worst, pair = 0.0, ""
        for i, ea in enumerate(envs):
            for eb in envs[i + 1:]:
                a = out[(out.env == ea) & (out.model == "XGBoost")]
                b = out[(out.env == eb) & (out.model == "XGBoost")]
                mm = a.merge(b, on=key, suffixes=("_a", "_b"))
                d = (mm[f"{c}_a"] - mm[f"{c}_b"]).abs().max()
                if d > worst:
                    worst, pair = d, f"{ea} vs {eb}"
        print(f"  {c:<8} {worst:.4f}  ({pair})")


if __name__ == "__main__":
    main(sys.argv[1:])
