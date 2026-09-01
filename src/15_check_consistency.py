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

    for b in bad:
        print(f"FAIL  {b}")
    if bad:
        sys.exit(f"{len(bad)} consistency failure(s)")
    print(f"consistency: {len(checks)} headline numbers agree across README, "
          f"METHODS.md and facts.json")


if __name__ == "__main__":
    main()
