"""Fig. 1D - EdU incorporation, nT vs rBM: two-tailed unpaired Student's t-test (equal variances).

Input: CSV with columns replicate, nT, rBM (one value per quantified field; 12 per group).
Welch's t-test (unequal variances) is printed for reference only.
Usage: python fig1D_edu_ttest.py <edu_quantification.csv> [<out.txt>]
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

d = pd.read_csv(sys.argv[1])
nt, rbm = d["nT"].dropna().to_numpy(float), d["rBM"].dropna().to_numpy(float)
res = stats.ttest_ind(nt, rbm, equal_var=True)
ci = res.confidence_interval(0.95)
welch = stats.ttest_ind(nt, rbm, equal_var=False)
lines = []
for name, x in (("nT", nt), ("rBM", rbm)):
    lines.append(f"{name}: n = {len(x)}, mean = {x.mean():.4f}, SD = {x.std(ddof=1):.4f}, "
                 f"SEM = {stats.sem(x):.4f}")
lines += [f"Student's t-test (two-tailed, equal variances): t = {res.statistic:.4f}, "
          f"df = {len(nt) + len(rbm) - 2}, P = {res.pvalue:.4f}",
          f"Difference of means (nT - rBM) = {nt.mean() - rbm.mean():.4f}, "
          f"95% CI [{ci.low:.4f}, {ci.high:.4f}]",
          f"(Welch's t-test, for reference: P = {welch.pvalue:.4f})"]
print("\n".join(lines))
if len(sys.argv) > 2:
    open(sys.argv[2], "w").write("\n".join(lines) + "\n")
