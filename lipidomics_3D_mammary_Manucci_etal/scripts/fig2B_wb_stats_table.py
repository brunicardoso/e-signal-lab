"""
Build a summary table of the parametric repeated-measures statistics.

Parses the same `statistical_analysis_summary_*` files that the plots use
(via fig2B_wb_barplots.REGISTRY / stats_file), so the table is guaranteed consistent with the
figures. For every marker it reports the RM-ANOVA omnibus p-value and all 10
pairwise comparisons (FDR-corrected p, stars, significance flag), with the two
within-timepoint rBM-/rBM+ comparisons (the ones annotated on the plots) marked.

Outputs:
    <WB_OUT_DIR>/stats_summary_table.csv   (machine-readable, one row per marker x pair)
and prints a Markdown table to stdout.

Usage:
    WB_OUT_DIR=../outputs/fig2B python fig2B_wb_stats_table.py
"""

import csv
import itertools
import os

from fig2B_wb_barplots import (
    OUT_DIR,
    REGISTRY,
    stats_file,
    WITHIN_TP_COMPARISONS,
    parse_pairwise,
    parse_main_test_p,
    p_to_stars,
)

CONDITION_ORDER = ["T0", "nT48", "rBM48", "nT72", "rBM72"]
PAIR_ORDER = list(itertools.combinations(CONDITION_ORDER, 2))
WITHIN_TP = {frozenset(c) for c in WITHIN_TP_COMPARISONS}
OUT_CSV = os.path.join(OUT_DIR, "stats_summary_table.csv")

HEADER = ["Marker", "Comparison", "Within-timepoint -/+",
          "p (FDR)", "Stars", "Significant", "RM-ANOVA p"]


def collect_rows():
    rows = []
    for key in REGISTRY:
        sf = stats_file(key)
        pairwise = parse_pairwise(sf)
        anova_p = parse_main_test_p(sf)
        for a, b in PAIR_ORDER:
            p, is_sig = pairwise[frozenset((a, b))]
            stars = p_to_stars(p) if is_sig else "ns"
            rows.append({
                "Marker": key,
                "Comparison": f"{a} vs {b}",
                "Within-timepoint -/+": "yes" if frozenset((a, b)) in WITHIN_TP else "",
                "p (FDR)": f"{p:.4f}",
                "Stars": stars,
                "Significant": "yes" if is_sig else "no",
                "RM-ANOVA p": f"{anova_p:.4f}" if anova_p is not None else "NA",
            })
    return rows


def write_csv(rows):
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)
    print(f"saved {OUT_CSV}\n")


def print_markdown(rows):
    print("| " + " | ".join(HEADER) + " |")
    print("|" + "|".join("---" for _ in HEADER) + "|")
    last = None
    for r in rows:
        marker = r["Marker"] if r["Marker"] != last else ""
        anova = r["RM-ANOVA p"] if r["Marker"] != last else ""
        last = r["Marker"]
        print(f"| {marker} | {r['Comparison']} | {r['Within-timepoint -/+']} "
              f"| {r['p (FDR)']} | {r['Stars']} | {r['Significant']} | {anova} |")


def main():
    rows = collect_rows()
    write_csv(rows)
    print_markdown(rows)


if __name__ == "__main__":
    main()
