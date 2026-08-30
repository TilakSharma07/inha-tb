"""Step 16 - programmatic quality gate on every figure.

Text collisions and off-canvas labels are the two figure defects that survive a casual look
at a 300-dpi PNG: at that size a legend sitting on the data or a clipped tick label is a few
pixels and easy to miss. Rendering each figure and measuring the text bounding boxes catches
both mechanically. Exits non-zero on any defect, so `run_all.sh` fails rather than shipping
a figure with a label on top of the data.

Limitation worth stating: this checks text against text. A label overlapping a *bar* or a
*point* has no text bounding box to collide with and is invisible here - those still need an
eye on the rendered figure.
"""
import itertools, os, runpy, sys
import matplotlib
matplotlib.use("Agg")

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.join(D, "src")
MIN_OVERLAP_PX = 4.0


def visible_ticks(ax, axis):
    """Tick labels matplotlib will actually draw - those inside the axis limits.

    Out-of-range tick Text objects survive a set_xlim/set_ylim call and still report a
    window extent, so a checker that reads them reports collisions for labels no reader
    can see.
    """
    lo, hi = sorted(ax.get_xlim() if axis == "x" else ax.get_ylim())
    labs = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    locs = ax.get_xticks() if axis == "x" else ax.get_yticks()
    return [t for t, v in zip(labs, locs) if lo - 1e-9 <= v <= hi + 1e-9]


def text_items(fig):
    out = []
    for ax in fig.axes:
        cands = (ax.texts + [ax.title, ax.xaxis.label, ax.yaxis.label]
                 + visible_ticks(ax, "x") + visible_ticks(ax, "y"))
        out += [(ax, t) for t in cands if t.get_text().strip()]
        lg = ax.get_legend()
        if lg:
            out += [(ax, t) for t in lg.get_texts()]
            if lg.get_title().get_text().strip():
                out.append((ax, lg.get_title()))
    return out


def check(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()
    T = text_items(fig)
    bad = []
    for (a1, t1), (a2, t2) in itertools.combinations(T, 2):
        if a1 is not a2:
            continue
        b1, b2 = t1.get_window_extent(r), t2.get_window_extent(r)
        if not b1.overlaps(b2):
            continue
        area = (max(0, min(b1.x1, b2.x1) - max(b1.x0, b2.x0))
                * max(0, min(b1.y1, b2.y1) - max(b1.y0, b2.y0)))
        if area > MIN_OVERLAP_PX:
            bad.append(f"{name}: text overlap {t1.get_text()[:24]!r} x {t2.get_text()[:24]!r}")
    for _, t in T:
        b = t.get_window_extent(r)
        if b.x0 < -2 or b.y0 < -2 or b.x1 > W + 2 or b.y1 > H + 2:
            bad.append(f"{name}: off-canvas {t.get_text()[:28]!r}")
    return bad


def main():
    figs, bad = [], []
    ns = runpy.run_path(f"{SRC}/10_figures.py", run_name="not_main")
    figs += [("figure1", ns["figure1"]()), ("figure2", ns["figure2"]())]
    j = f"{D}/data/docking_joined.csv"
    if os.path.exists(j):
        ns3 = runpy.run_path(f"{SRC}/13_figure_docking.py", run_name="not_main")
        figs.append(("figure3", ns3["main"]()))
    else:
        print("figure3 skipped (data/docking_joined.csv absent)")

    for name, fig in figs:
        errs = check(fig, name)
        print(f"{name}: {len(errs)} defect(s)")
        bad += errs
    for b in bad:
        print(f"FAIL  {b}")
    if bad:
        sys.exit(f"{len(bad)} figure defect(s)")
    print(f"figures: {len(figs)} checked, no text collisions or clipped labels")


if __name__ == "__main__":
    main()
