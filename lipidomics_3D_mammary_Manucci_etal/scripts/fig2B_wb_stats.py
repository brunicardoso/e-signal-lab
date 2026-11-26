"""
Parametric repeated-measures statistical analysis of Western-blot z-scores.

For every marker (ACC, pACC, pACC/ACC, CPT1A, CPT2, LamB1) this analyses the
`Z_score` column with:

    Main test : Repeated measures ANOVA (statsmodels AnovaRM)
    Post-hoc  : paired t-tests (all condition pairs), with
                Benjamini-Hochberg FDR correction across the pairwise family

Why parametric: with n=4 replicates, rank-based tests (Wilcoxon/Friedman+Nemenyi)
are mathematically unable to detect the consistent low-variance effects in this
data (Wilcoxon's smallest two-sided p is 0.125). The z-scores pass Shapiro-Wilk
normality for essentially all groups and the design is within-subjects, so the
parametric paired approach is both valid and adequately powered.

The report layout mirrors the original parametric summaries
(statistical_analysis_summary_20260519_160347_ACC.txt) with an added
"Median [IQR]" descriptive line.

Output: one file per marker in <out-dir>
    statistical_analysis_summary_parametric_<timestamp>_<marker>.txt

Usage:
    python fig2B_wb_stats.py --data-dir ../data/western_blot --out-dir ../outputs/fig2B
"""

import argparse
import os
from datetime import datetime
import itertools
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
VALUE_COL = "Z_score"
ALPHA = 0.05
CONDITION_ORDER = ["T0", "nT48", "rBM48", "nT72", "rBM72"]

FILES = [
    ("ACC", "ACC_CellBioStats_version.xlsx"),
    ("pACC", "pACC_CellBioStats_version.xlsx"),
    ("pACC_per_ACC", "pACC_per_ACC_CellBioStats_version.xlsx"),
    ("CPT1A", "CPT1A_CellBioStats_version.xlsx"),
    ("CPT2", "CPT2_CellBioStats_version.xlsx"),
    ("LamB1", "LamB1_CellBioStats_version.xlsx"),
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def build_matrix(df):
    """(blocks x groups) DataFrame: rows=Replicate, cols=CONDITION_ORDER."""
    mat = df.pivot(index="Replicate", columns="Condition", values=VALUE_COL)
    missing = [c for c in CONDITION_ORDER if c not in mat.columns]
    if missing:
        raise ValueError(f"missing conditions {missing}")
    mat = mat[CONDITION_ORDER]
    if mat.isna().any().any():
        raise ValueError("incomplete blocks (NaN present) - RM-ANOVA needs a "
                         "balanced design")
    return mat


def descriptives(values):
    vals = np.asarray(values, dtype=float)
    n = len(vals)
    mean = vals.mean()
    sem = vals.std(ddof=1) / np.sqrt(n) if n > 1 else 0.0
    median = np.median(vals)
    q1, q3 = np.percentile(vals, [25, 75])
    return mean, sem, median, q1, q3, n


def rm_anova_p(df):
    """Repeated-measures ANOVA p-value (Pr > F) for Condition."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = AnovaRM(df, depvar=VALUE_COL, subject="Replicate",
                      within=["Condition"]).fit()
    return float(res.anova_table["Pr > F"].iloc[0])


def pairwise_fdr(mat):
    """Paired t-tests for all condition pairs with BH-FDR correction.

    Returns ordered list of (cond_i, cond_j, p_corrected).
    """
    pairs = list(itertools.combinations(CONDITION_ORDER, 2))
    raw = [stats.ttest_rel(mat[a], mat[b]).pvalue for a, b in pairs]
    corrected = multipletests(raw, method="fdr_bh")[1]
    return [(a, b, float(p)) for (a, b), p in zip(pairs, corrected)]


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def build_report(df, mat):
    lines = []
    lines.append("Statistical Analysis Summary")
    lines.append("=" * 50)
    lines.append("")
    lines.append("Note: Parametric repeated-measures test (RM-ANOVA) was used.")
    lines.append("")

    # descriptives
    lines.append("Descriptive Statistics (based on replicate means):")
    lines.append("-" * 50)
    for cond in CONDITION_ORDER:
        mean, sem, median, q1, q3, n = descriptives(mat[cond].to_numpy(float))
        lines.append(f"Group {cond}:")
        lines.append(f"  Mean +/- SEM: {mean:.2f} +/- {sem:.2f}")
        lines.append(f"  Median [IQR]: {median:.2f} [{q1:.2f}, {q3:.2f}]")
        lines.append(f"  N (replicates) = {n}")
        lines.append(f"  N (total observations) = {n}")
    lines.append("")

    # normality
    lines.append("Normality Test (Shapiro-Wilk):")
    lines.append("-" * 50)
    for cond in CONDITION_ORDER:
        p = stats.shapiro(mat[cond].to_numpy(float)).pvalue
        lines.append(f"Group {cond}:")
        lines.append(f"  p-value: {p:.4f} ({'Normal' if p > ALPHA else 'Not normal'})")
    lines.append("")

    # homoscedasticity
    lines.append("Homoscedasticity Test (Levene):")
    lines.append("-" * 50)
    groups = [mat[c].to_numpy(float) for c in CONDITION_ORDER]
    lev_p = stats.levene(*groups).pvalue
    lines.append(f"p-value: {lev_p:.4f}")
    lines.append("(Homoscedastic)" if lev_p > ALPHA else "(Heteroscedastic)")
    lines.append("")

    # main test
    lines.append("Main Test (Repeated measures ANOVA (AnovaRM)):")
    lines.append("-" * 50)
    lines.append(f"p-value: {rm_anova_p(df):.4f}")
    lines.append("")

    # pairwise
    lines.append("Pairwise Comparisons:")
    lines.append("-" * 50)
    for a, b, p in pairwise_fdr(mat):
        verdict = "Significant" if p < ALPHA else "Not significant"
        lines.append(f"{a} vs {b}: p = {p:.4f} ({verdict}) "
                     f"[Paired t-test (Benjamini-Hochberg FDR)]")

    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=".", help="folder with the *_CellBioStats_version.xlsx files")
    ap.add_argument("--out-dir", default=".", help="folder for the summary .txt files")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for key, xlsx in FILES:
        df = pd.read_excel(os.path.join(args.data_dir, xlsx))
        mat = build_matrix(df)
        report = build_report(df, mat)
        out = os.path.join(args.out_dir,
                           f"statistical_analysis_summary_parametric_{timestamp}_{key}.txt")
        with open(out, "w") as fh:
            fh.write(report)
        n_sig = sum(p < ALPHA for *_, p in pairwise_fdr(mat))
        print(f"{key}: RM-ANOVA p = {rm_anova_p(df):.4f}  |  "
              f"{n_sig} FDR-significant pairs  ->  {out}")


if __name__ == "__main__":
    main()
