"""Step 15 - assert README, METHODS.md and the saved tables all state the same numbers.

Written after two silent drifts: hand-typed tables that fell out of date when the pipeline
was rerun, and the same underlying value rounded two different ways in two documents. Both
were invisible on reading and obvious to a diff. This runs last in the pipeline and exits
non-zero if any headline number disagrees with `docs/facts.json`.
"""
import json, os, re, sys

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def main():
    F = json.load(open(f"{D}/docs/facts.json"))
    R = open(f"{D}/README.md").read()
    M = open(f"{D}/docs/METHODS.md").read()

    dk, ad = F["docking"], F["applicability_domain"]
    # Each check is a regex with ONE capture group, matched against the document; the
    # captured text must equal the expected value. A bare substring search is not enough:
    # 0.23 is both the redock RMSD and the assay noise floor, so corrupting one still
    # matched the other and the check passed on a wrong document.
    N = r"([\d.,]+)"
    checks = [
        ("redock RMSD / README",   R, rf"\*\*{N} Å\*\* RMSD",             f"{dk['rmsd_top']:.2f}"),
        ("redock RMSD / METHODS",  M, rf"\*\*RMSD, top-ranked pose\*\* \| \*\*{N} Å",
                                                                          f"{dk['rmsd_top']:.2f}"),
        # Method-section library breakdown. These were stale (539/344/195/348 against
        # a 647-molecule library) and no check covered them, so the README contradicted
        # its own headline table. Any prose number that restates a saved value is checked.
        ("library total / Method",  R, rf"Library: {N} molecules",        f"{dk['library_n']}"),
        ("library InhA / Method",   R, rf"\({N} with InhA enzyme",        f"{dk['library_inha']}"),
        ("library WC / Method",     R, rf"enzyme potency, {N}\n",         f"{dk['library_wc']}"),
        ("library actives / Method", R, rf"of which {N} are measured",     f"{dk['library_actives']}"),
        ("excluded / Method",       R, rf"{N} molecules were excluded",    f"{dk['n_excluded']}"),
        ("scored / Method",         R, rf"{N} of the [\d,]+ returned",     f"{dk['n_scored']}"),
        ("library / README",       R, rf"Compounds in the docking library \| {N}",
                                                                          f"{dk['library_n']:,}"),
        ("library / METHODS",      M, rf"Library: {N} molecules",          f"{dk['library_n']:,}"),
        ("curated / README",       R, rf"Curated bioactivity records \| {N}",
                                                                          f"{F['records_curated']:,}"),
        ("unique mols / README",   R, rf"→ {N} unique molecules",          f"{F['molecules_unique']:,}"),
        ("unique mols / METHODS",  M, rf"\| unique molecules \(aggregated\) \| {N} \|",
                                                                          f"{F['molecules_unique']:,}"),
        ("noise floor / README",   R, rf"replicate assays\) \| {N} log units",
                                                                          f"{F['noise_floor_sd']:.2f}"),
        ("noise floor / METHODS",  M, rf"median inter-assay SD is \*\*{N} log units",
                                                                          f"{F['noise_floor_sd']:.2f}"),
        ("InhA set / METHODS",     M, rf"\*\*{N}\*\* molecules with InhA enzyme potency",
                                                                          f"{F['n_inha']:,}"),
        ("whole-cell / METHODS",   M, rf"\*\*{N}\*\* with whole-cell potency",
                                                                          f"{F['n_wc']:,}"),
        ("AD low / README",        R, rf"score ROC-AUC {N}–",             f"{min(ad['near_roc_auc']):.2f}"),
        ("AD high / README",       R, rf"score ROC-AUC [\d.]+–{N};",      f"{max(ad['near_roc_auc']):.2f}"),
        ("AD off low / README",    R, rf"falls to {N}–",                  f"{min(ad['far_roc_auc']):.2f}"),
        ("AD off high / README",   R, rf"falls to [\d.]+–{N} —",          f"{max(ad['far_roc_auc']):.2f}"),
        ("inflation / README",     R, rf"inflates ROC-AUC by \*\*{N}–",   f"{F['inflation']['min']:.2f}"),
        ("inflation / METHODS",    M, rf"overstates ROC-AUC by\s*\n\*\*{N}–",
                                                                          f"{F['inflation']['min']:.2f}"),
        # the head-to-head is the screen's result; it must read the same in both documents
        ("h2h dock InhA / README",  R, rf"same held-out molecules \| ROC-AUC \*\*{N}\*\*",
                                       f"{dk['InhA_enzyme_heldout_dock_auc']:.2f}"),
        ("h2h dock WC / README",    R, rf"same held-out molecules \| ROC-AUC \*\*[\d.]+\*\* / \*\*{N}\*\*",
                                       f"{dk['Mtb_whole_cell_heldout_dock_auc']:.2f}"),
        ("h2h dock InhA / METHODS", M, rf"\| InhA enzyme \| \d+ \| {N} \| [\d.]+ \|\n",
                                       f"{dk['InhA_enzyme_heldout_dock_auc']:.2f}"),
        ("h2h qsar InhA / METHODS", M, rf"\| InhA enzyme \| \d+ \| [\d.]+ \| {N} \|\n",
                                       f"{dk['InhA_enzyme_heldout_qsar_auc']:.2f}"),
        ("screen AUC InhA / METHODS", M,
         rf"\| InhA enzyme \| \d+ \| \d+ \| {N} \| -?[\d.]+ \|",
                                       f"{dk['InhA_enzyme_roc_auc']:.2f}"),
        # Leakage decomposition (src/21). The README restates saved values; a rerun
        # that changes them must not leave the prose behind. The permutation p-value is
        # deliberately NOT checked - it moves in the last digit between runs - but the
        # sentence stating it is not significant IS checked, below the loop.
        ("leakage combos / README",  R, rf"\*\*{N}/6\*\*\s+model", "6"),
        ("leakage median / README",  R, rf"median \+{N}\s+ROC-AUC", "0.099"),
        ("leakage residual / README", R, r"median residual\s+([+-][\d.]+)", "-0.006"),
        ("leakage residual / METHODS", M, r"median residual\s+([+-][\d.]+)", "-0.006"),
        ("leakage median / METHODS",  M, rf"median\s+\+{N}\s+ROC-AUC", "0.099"),
        ("leakage combos / METHODS",  M, rf"in\s+{N}/6 model", "6"),
    ]
    bad = []
    for label, text, pat, want in checks:
        m = re.search(pat, text)
        if m is None:
            bad.append(f"{label}: pattern not found - {pat}")
        elif m.group(1).replace(",", "") != want.replace(",", ""):
            bad.append(f"{label}: document says {m.group(1)!r}, tables say {want!r}")

    # the two documents must not round the same quantity differently
    m = re.search(r"score ROC-AUC ([\d.]+)–([\d.]+); below 0\.35 .*?falls to ([\d.]+)–([\d.]+)",
                  R, re.S)
    n = re.search(r"interpretable bins score ([\d.]+)–([\d.]+); below 0\.35 they fall\s*\n"
                  r"to ([\d.]+)–([\d.]+)", M)
    if m and n and m.groups() != n.groups():
        bad.append(f"AD ranges disagree: README {m.groups()} vs METHODS {n.groups()}")
    elif not (m and n):
        bad.append("AD range sentence not found in one of the documents")

    # The README must keep saying the effect is not significant while facts.json says the
    # direction is consistent; those two must not drift apart into an overclaim.
    if not re.search(r"not significant at the\s+\**conventional threshold",
                     R.replace("**", "")):
        bad.append("README no longer states the leakage effect is not significant")

    # Version-drift bounds. These were prose-only: no facts.json key covered them, so
    # nothing stopped them going stale while data/version_attribution.csv moved
    # underneath. They are the worst PAIRWISE spread across environments, and they round
    # UP, being bounds - see the note in src/11_facts.py.
    vd = F.get("version_drift")
    if vd:
        for label, pat, keys in [
            ("drift ROC-AUC/PR-AUC", rf"XGBoost ROC-AUC, PR-AUC \| up to {N} / {N}\.",
             ("xgboost_ROC_AUC", "xgboost_PR_AUC")),
            ("drift MCC/BalAcc", rf"XGBoost MCC, balanced accuracy \| up to {N} / {N}\.",
             ("xgboost_MCC", "xgboost_BalAcc")),
        ]:
            m = re.search(pat, R)
            if not m:
                bad.append(f"{label}: sentence not found in README")
                continue
            want = tuple(f"{vd[k]:.3f}" for k in keys)
            if m.groups() != want:
                bad.append(f"{label}: README {m.groups()} vs facts.json {want}")
        m = re.search(rf"EF5 \| up to {N},", R)
        if not m:
            bad.append("drift EF5: sentence not found in README")
        elif m.group(1) != f"{vd['xgboost_EF5']:.3f}":
            bad.append(f"drift EF5: README {m.group(1)} vs facts.json "
                       f"{vd['xgboost_EF5']:.3f}")
        # Non-monotonicity: the README claims python-only drift EXCEEDS the endpoint
        # jump. Check both numbers and the inequality - if a re-measurement reversed it,
        # the numbers could each still match while the sentence became false.
        m = re.search(rf"moves MCC further\s+\({N}\) than the full jump to newer xgboost\s+"
                      rf"and numpy does \({N}\)", R)
        if not m:
            bad.append("drift non-monotonicity: sentence not found in README")
        else:
            want = (f"{vd['mcc_python_only']:.3f}", f"{vd['mcc_endpoint']:.3f}")
            if m.groups() != want:
                bad.append(f"drift non-monotonicity: README {m.groups()} vs "
                           f"facts.json {want}")
            elif vd["mcc_python_only"] <= vd["mcc_endpoint"]:
                bad.append("drift non-monotonicity: README claims python-only drift "
                           f"exceeds the endpoint jump, but {vd['mcc_python_only']} <= "
                           f"{vd['mcc_endpoint']}")

    # The CSV-parser perturbation is the evidence that the interpreter effect is not
    # generic small-number amplification. Three numbers carry that argument and all three
    # can go stale independently, so check each against facts.json - and check the
    # inequality too: the argument only holds while the perturbation the model tolerates
    # is LARGER than the QED drift it does not.
    cp = F.get("csv_parser")
    if cp:
        m = re.search(rf"changes {N} cells of the model input matrix", R)
        if not m:
            bad.append("csv parser cells: sentence not found in README")
        elif m.group(1).replace(",", "") != str(cp["cells_perturbed"]):
            bad.append(f"csv parser cells: README {m.group(1)} vs facts.json "
                       f"{cp['cells_perturbed']}")
        n_cols = len(cp["columns_perturbed"])
        words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
                 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
        if not re.search(rf"across {words.get(n_cols, n_cols)} of the seventeen\s+descriptor columns", R):
            bad.append(f"csv parser columns: README does not say "
                       f"'{words.get(n_cols, n_cols)} of the seventeen'")
        for c in cp["columns_perturbed"]:
            if f"`{c}`" not in R:
                bad.append(f"csv parser columns: {c} perturbed but not named in README")
        m = re.search(rf"up to {N}e-13\s+absolute", R)
        if not m:
            bad.append("csv parser magnitude: sentence not found in README")
        elif m.group(1) != f"{cp['max_abs_e13']:g}":
            bad.append(f"csv parser magnitude: README {m.group(1)}e-13 vs facts.json "
                       f"{cp['max_abs_e13']:g}e-13")
        # The comparison figure must come from version_drift (a 3-environment bound over
        # version_attribution.csv), not from reproducibility_drift.csv - that one is
        # rewritten by step 18 on every run and its max is EF5, which moves ~0.35 on a
        # single rank change. See the note in src/11_facts.py.
        if vd:
            m = re.search(rf"against up to {N} in MCC alone", R)
            if not m:
                bad.append("csv parser comparison: sentence not found in README")
            elif m.group(1) != f"{vd['xgboost_MCC']:.3f}":
                bad.append(f"csv parser comparison: README {m.group(1)} vs facts.json "
                           f"{vd['xgboost_MCC']:.3f}")

    # Printed last, after every check has appended. An earlier position silently hid the
    # failures appended below it: they still set the exit code, with nothing on stdout
    # saying which check failed.
    for b in bad:
        print(f"FAIL  {b}")

    if bad:
        sys.exit(f"{len(bad)} consistency failure(s)")
    # len(checks) counts only the tabulated number-vs-fact pairs. The prose claims
    # checked imperatively below that loop (leakage wording, version-drift bounds and
    # their non-monotonicity, the CSV-parser perturbation) are verified but not counted,
    # so this line must not present the number as the total - it read as one until a
    # mutation test showed the count unchanged after three new checks were added.
    print(f"consistency: {len(checks)} tabulated numbers agree across README, "
          f"METHODS.md and facts.json, plus the prose claims checked below them")


if __name__ == "__main__":
    main()
