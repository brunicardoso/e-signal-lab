"""
Bar plots of Western-blot quantification, styled like the reference figure.

Works with any "*_CellBioStats_version.xlsx" file that has columns
    Condition | Measument | Replicate | Z_score | Min_Max_Norm
and the matching "statistical_analysis_summary_*.txt" file.

For each input it produces two plots (one per value column), PNG/SVG/PDF:
    <prefix>_Zscore_barplot.*
    <prefix>_MinMaxNorm_barplot.*

Condition order / colours follow the reference image:
    T0 (black) | 48 h: nT48 (teal, rBM-), rBM48 (green, rBM+) | 72 h: nT72, rBM72

Significance stars are read automatically from the stats file's
"Pairwise Comparisons" section (only significant pairs are drawn,
GraphPad thresholds).

Usage
-----
Plot every file in the built-in registry (default):
    WB_DATA_DIR=../data/western_blot WB_OUT_DIR=../outputs/fig2B python fig2B_wb_barplots.py

Plot a single custom file:
    python fig2B_wb_barplots.py \
        --xlsx MyProtein_CellBioStats_version.xlsx \
        --stats statistical_analysis_summary_xxx.txt \
        --label "CPT1A" --prefix CPT1A
"""

import argparse
import glob
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# Layout shared by all plots (edit here to restyle every figure at once)
# --------------------------------------------------------------------------- #
# key, x-position, colour, rBM label, time label
CONDITIONS = [
    ("T0",    0.0,  "#000000", "-", "0"),
    ("nT48",  1.10, "#2c6e7f", "-", "48"),
    ("rBM48", 1.95, "#82c341", "+", "48"),
    ("nT72",  3.05, "#2c6e7f", "-", "72"),
    ("rBM72", 3.90, "#82c341", "+", "72"),
]
BAR_WIDTH = 0.75

# Value columns to plot -> (column name, prefix suffix, y-axis sub-label)
VALUE_COLUMNS = [
    ("Z_score", "Zscore", "(z-score)"),
    ("Min_Max_Norm", "MinMaxNorm", "(min-max norm)"),
]


# --------------------------------------------------------------------------- #
# Stats parsing
# --------------------------------------------------------------------------- #
# e.g. "T0 vs rBM48: p = 0.0008 (Significant) [Paired t-test (Bonferroni)]"
_PAIR_RE = re.compile(
    r"^\s*(\S+)\s+vs\s+(\S+)\s*:\s*p\s*=\s*([0-9.eE+-]+)\s*\((Significant|Not significant)\)"
)


# Main-test "p-value:" line (used by the summary-table script)
_MAIN_P_RE = re.compile(r"^p-value:\s*([0-9.eE+-]+)")

# Comparisons annotated on the plots: within-timepoint rBM- (nT) vs rBM+ (rBM).
WITHIN_TP_COMPARISONS = [("nT48", "rBM48"), ("nT72", "rBM72")]


def parse_pairwise(stats_path):
    """Return {frozenset({a, b}): (p_value, is_significant)} for ALL pairs."""
    out = {}
    with open(stats_path) as fh:
        in_section = False
        for line in fh:
            if line.startswith("Pairwise Comparisons"):
                in_section = True
                continue
            if in_section:
                m = _PAIR_RE.match(line)
                if m:
                    a, b, p, flag = m.groups()
                    out[frozenset((a, b))] = (float(p), flag == "Significant")
    return out


def parse_main_test_p(stats_path):
    """Return the Main Test (omnibus) p-value."""
    with open(stats_path) as fh:
        in_section = False
        for line in fh:
            if line.startswith("Main Test"):
                in_section = True
                continue
            if in_section:
                m = _MAIN_P_RE.match(line.strip())
                if m:
                    return float(m.group(1))
    return None


def p_to_stars(p):
    """GraphPad-style significance stars."""
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def _group_stats(df, value_col):
    """condition -> (values, mean, sem)."""
    out = {}
    for key, *_ in CONDITIONS:
        vals = df.loc[df["Condition"] == key, value_col].to_numpy(dtype=float)
        mean = vals.mean()
        sem = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0
        out[key] = (vals, mean, sem)
    return out


def _make_plot(df, value_col, ylabel_sub, protein_label, out_stem, pairwise):
    stats = _group_stats(df, value_col)
    xpos = {c[0]: c[1] for c in CONDITIONS}
    rng = np.random.default_rng(7)  # reproducible jitter

    fig, ax = plt.subplots(figsize=(5.2, 4.8))

    # bars, SEM error bars, individual replicate points
    for key, x, colour, _rbm, _time in CONDITIONS:
        vals, mean, sem = stats[key]
        ax.bar(x, mean, width=BAR_WIDTH, color=colour, zorder=2)
        ax.errorbar(x, mean, yerr=sem, fmt="none", ecolor="black",
                    elinewidth=1.4, capsize=5, capthick=1.4, zorder=4)
        jitter = (rng.random(len(vals)) - 0.5) * BAR_WIDTH * 0.45
        ax.scatter(np.full(len(vals), x) + jitter, vals,
                   facecolors="white", edgecolors="black", linewidths=1.1,
                   s=34, zorder=5)

    # data extent for placing brackets
    tops, bots = [], []
    for key in xpos:
        vals, mean, sem = stats[key]
        tops.append(max(vals.max(), mean + sem))
        bots.append(min(vals.min(), mean - sem))
    data_top, data_bot = max(tops), min(bots)
    yspan = data_top - data_bot if data_top > data_bot else 1.0

    # within-timepoint rBM- vs rBM+ brackets (always drawn; 'ns' if not sig).
    # each sits just above the taller of its two bars (groups don't overlap in x)
    bar_h = yspan * 0.02
    for a, b in WITHIN_TP_COMPARISONS:
        if a not in xpos or b not in xpos or frozenset((a, b)) not in pairwise:
            continue
        p, is_sig = pairwise[frozenset((a, b))]
        label = p_to_stars(p) if is_sig else "ns"
        x1, x2 = sorted((xpos[a], xpos[b]))
        local_top = max(max(stats[a][0].max(), stats[a][1] + stats[a][2]),
                        max(stats[b][0].max(), stats[b][1] + stats[b][2]))
        h = local_top + yspan * 0.08
        ax.plot([x1, x1, x2, x2], [h, h + bar_h, h + bar_h, h],
                color="black", linewidth=1.2, zorder=6)
        ax.text((x1 + x2) / 2, h + bar_h, label,
                ha="center", va="bottom", fontsize=15, zorder=6)

    # axes styling
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(1.6)
    ax.tick_params(axis="y", width=1.6, length=6, labelsize=14)
    ax.tick_params(axis="x", length=0)
    ax.axhline(0, color="black", linewidth=1.6, zorder=1)

    # y label: bold protein name + smaller score sub-label
    ax.set_ylabel("")
    ax.text(-0.34, 0.5, protein_label, transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=30, fontweight="bold")
    ax.text(-0.22, 0.5, ylabel_sub, transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=15)

    ax.set_xticks([])
    ax.set_xlim(-0.7, CONDITIONS[-1][1] + 0.7)

    # annotation rows below the axis (rBM / Time)
    ymin, ymax = ax.get_ylim()
    row1 = ymin - 0.10 * (ymax - ymin)
    row2 = ymin - 0.20 * (ymax - ymin)
    ax.set_ylim(ymin, ymax)

    for key, x, colour, rbm, time in CONDITIONS:
        ax.text(x, row1, rbm, ha="center", va="center", fontsize=18, clip_on=False)
    ax.text(-0.65, row1, "rBM", ha="right", va="center", fontsize=15, clip_on=False)
    ax.text(-0.65, row2, "Time (h)", ha="right", va="center", fontsize=15, clip_on=False)

    groups = {}
    for key, x, colour, rbm, time in CONDITIONS:
        groups.setdefault(time, []).append(x)
    for time, xs in groups.items():
        ax.text(np.mean(xs), row2, time, ha="center", va="center",
                fontsize=18, clip_on=False)
        if len(xs) > 1:
            yline = (row1 + row2) / 2
            ax.plot([min(xs) - 0.1, max(xs) + 0.1], [yline, yline],
                    color="black", linewidth=1.3, clip_on=False)

    plt.tight_layout()
    fig.subplots_adjust(left=0.30)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_stem}.png / .svg / .pdf")


def plot_wb(xlsx, stats_file, protein_label, prefix):
    """Generate both score plots for one protein."""
    df = pd.read_excel(xlsx)
    pairwise = parse_pairwise(stats_file)
    shown = []
    for a, b in WITHIN_TP_COMPARISONS:
        if frozenset((a, b)) in pairwise:
            p, is_sig = pairwise[frozenset((a, b))]
            shown.append(f"{a} vs {b}: {p_to_stars(p) if is_sig else 'ns'}")
    print(f"{prefix}: {' | '.join(shown)}")
    for value_col, suffix, ylabel_sub in VALUE_COLUMNS:
        _make_plot(df, value_col, ylabel_sub, protein_label,
                   f"{prefix}_{suffix}_barplot", pairwise)


# --------------------------------------------------------------------------- #
# Registry of known files: prefix -> (xlsx, stats_file, y-axis protein label)
# --------------------------------------------------------------------------- #
# Input/output folders: environment variables WB_DATA_DIR / WB_OUT_DIR (default: current folder).
# The stats file of each marker is the newest
# statistical_analysis_summary_parametric_<timestamp>_<marker>.txt in WB_OUT_DIR
# (written by fig2B_wb_stats.py).
DATA_DIR = os.environ.get("WB_DATA_DIR", ".")
OUT_DIR = os.environ.get("WB_OUT_DIR", ".")
REGISTRY = {
    "ACC": ("ACC_CellBioStats_version.xlsx", "ACC"),
    "pACC": ("pACC_CellBioStats_version.xlsx", "p-ACC"),
    "pACC_per_ACC": ("pACC_per_ACC_CellBioStats_version.xlsx", "p-ACC/ACC"),
    "CPT1A": ("CPT1A_CellBioStats_version.xlsx", "CPT1A"),
    "CPT2": ("CPT2_CellBioStats_version.xlsx", "CPT2"),
    "LamB1": ("LamB1_CellBioStats_version.xlsx", "LamB1"),
}


def stats_file(key):
    """Newest parametric stats summary for `key` in OUT_DIR."""
    hits = sorted(glob.glob(os.path.join(
        OUT_DIR, f"statistical_analysis_summary_parametric_*_{key}.txt")))
    hits = [h for h in hits if h.endswith(f"_{key}.txt") and
            re.fullmatch(r"statistical_analysis_summary_parametric_\d{8}_\d{6}_" + re.escape(key) + r"\.txt",
                         os.path.basename(h))]
    if not hits:
        raise FileNotFoundError(f"no stats summary for {key} in {OUT_DIR}; run fig2B_wb_stats.py first")
    return hits[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", help="data .xlsx file")
    ap.add_argument("--stats", help="statistical_analysis_summary_*.txt file")
    ap.add_argument("--label", help="bold y-axis protein label (e.g. CPT1A)")
    ap.add_argument("--prefix", help="output filename prefix")
    ap.add_argument("--only", nargs="*", choices=list(REGISTRY),
                    help="restrict registry run to these keys")
    args = ap.parse_args()

    if args.xlsx:  # single custom file
        if not (args.stats and args.label and args.prefix):
            ap.error("--xlsx requires --stats, --label and --prefix")
        plot_wb(args.xlsx, args.stats, args.label, args.prefix)
        return

    keys = args.only or list(REGISTRY)
    os.makedirs(OUT_DIR, exist_ok=True)
    for key in keys:
        xlsx, label = REGISTRY[key]
        plot_wb(os.path.join(DATA_DIR, xlsx), stats_file(key), label, os.path.join(OUT_DIR, key))


if __name__ == "__main__":
    main()
