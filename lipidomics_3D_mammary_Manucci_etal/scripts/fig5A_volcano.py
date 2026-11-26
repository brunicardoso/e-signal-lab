"""Fig. 5A - volcano plot, rBM vs nT, using log ratios.

Thresholds: FDR (BH) < 0.10 and |FC| > 1.5.
  log2FC = mean(log2 rBM) - mean(log2 nT)   (ratio of geometric means)
  p      = Student's t-test (equal var) on log10 values
Usage: python fig5A_volcano.py <NormDNA.csv> <out_prefix>
"""
import sys, csv
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV, OUT = sys.argv[1], sys.argv[2]
GROUP_NUM, GROUP_DEN, COND_LABEL = "lrECM", "nT", "rBM"
FC_THRESHOLD = 1.5
LOG2FC_CUTOFF = np.log2(FC_THRESHOLD)
FDR_THRESHOLD = 0.10
C_NS, C_UP, C_DN = "#CCCCCC", "#F0A30A", "#56B4E9"

rows = list(csv.reader(open(CSV, encoding="utf-8-sig")))
hdr = rows[0]
data = [dict(zip(hdr, r)) for r in rows[1:]]
species = [c for c in hdr if c not in ("Sample", "Group")]

def grp(g, sp):
    return np.array([float(r[sp]) for r in data if r["Group"] == g])

def bh(p):
    p = np.asarray(p); m = len(p); i = np.argsort(p)
    q = p[i] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    o = np.empty(m); o[i] = np.minimum(q, 1); return o

res, skipped = [], []
for sp in species:
    a, b = grp(GROUP_DEN, sp), grp(GROUP_NUM, sp)
    if (a <= 0).any() or (b <= 0).any():
        skipped.append(sp); continue
    la, lb = np.log10(a), np.log10(b)
    p = stats.ttest_ind(lb, la, equal_var=True).pvalue
    res.append(dict(Identifier=sp, nT_geomean=10 ** la.mean(), lrECM_geomean=10 ** lb.mean(),
                    log2FC=np.log2(10) * (lb.mean() - la.mean()), p_value=p))
fdr = bh([r["p_value"] for r in res])
for r, f in zip(res, fdr):
    r["FDR"] = f
    r["FC"] = 2 ** r["log2FC"]
    r["neg_log10_FDR"] = -np.log10(f)
    r["Direction"] = "Up" if r["log2FC"] > 0 else "Down"
    r["Significant"] = bool(f < FDR_THRESHOLD and abs(r["log2FC"]) > LOG2FC_CUTOFF)

print(f"Species tested: {len(res)}; skipped (zeros): {skipped}")
up = sorted([r for r in res if r["Significant"] and r["Direction"] == "Up"], key=lambda r: r["FDR"])
dn = sorted([r for r in res if r["Significant"] and r["Direction"] == "Down"], key=lambda r: r["FDR"])
ns = [r for r in res if not r["Significant"]]
print(f"Significant: {len(up) + len(dn)} (up={len(up)}, down={len(dn)})")

keys = ["Identifier", "nT_geomean", "lrECM_geomean", "FC", "log2FC", "p_value", "FDR",
        "neg_log10_FDR", "Direction", "Significant"]
with open(OUT + ".csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(res)

# --- plot (same style as volcano_plot_FDR010) ---
fig, ax = plt.subplots(figsize=(13, 8))
xy = lambda L: ([r["log2FC"] for r in L], [r["neg_log10_FDR"] for r in L])
ax.scatter(*xy(ns), c=C_NS, s=45, alpha=0.6, edgecolors="none", label="Not Significant")
ax.scatter(*xy(up), c=C_UP, s=70, edgecolors="none", label=f"Upregulated ({COND_LABEL} > {GROUP_DEN})")
ax.scatter(*xy(dn), c=C_DN, s=70, edgecolors="none", label=f"Downregulated ({COND_LABEL} < {GROUP_DEN})")
ax.axhline(-np.log10(FDR_THRESHOLD), ls="--", c="gray", lw=1)
ax.axvline(LOG2FC_CUTOFF, ls="--", c="gray", lw=1)
ax.axvline(-LOG2FC_CUTOFF, ls="--", c="gray", lw=1)

labelmap, n = [], 1
for direction, L in (("Up", up), ("Down", dn)):
    for r in L:
        ax.annotate(str(n), (r["log2FC"], r["neg_log10_FDR"]), textcoords="offset points",
                    xytext=(-5, 7), fontsize=10)
        labelmap.append((n, direction, r["Identifier"])); n += 1

ax.set_xlabel("log2(Fold Change)", fontsize=18, fontweight="bold")
ax.set_ylabel("-log10(FDR)", fontsize=18, fontweight="bold")
ax.legend(loc="lower left", frameon=True, fontsize=11)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", which="major", labelsize=17)

fig.subplots_adjust(right=0.74)
line_h, x0, y = 0.045, 1.03, 0.99
for title, col, d in (("Upregulated:", C_UP, "Up"), ("Downregulated:", C_DN, "Down")):
    ax.text(x0, y, title, transform=ax.transAxes, fontsize=14, fontweight="bold", color=col, va="top")
    y -= line_h * 1.3
    for num_, d_, name_ in [t for t in labelmap if t[1] == d]:
        ax.text(x0, y, f"{num_} = {name_.replace('1G-cer', '1G-Cer')}",   # display label only
                transform=ax.transAxes, fontsize=12, va="top")
        y -= line_h
    y -= line_h

fig.text(0.125, -0.02,
         "log2FC = log2 ratio of geometric means (mean of log2 rBM − mean of log2 nT). "
         "Student's t-test on log10 data; Benjamini–Hochberg FDR across all species.\n"
         f"Significant: FDR < {FDR_THRESHOLD:g} and |FC| > {FC_THRESHOLD:g}; n = 4 per group.",
         fontsize=10, color="0.3", va="top")
plt.savefig(OUT + ".png", dpi=300, bbox_inches="tight")
plt.savefig(OUT + ".svg", bbox_inches="tight")
for t in labelmap:
    r = next(x for x in res if x["Identifier"] == t[2])
    print(t, f"log2FC={r['log2FC']:+.2f} FDR={r['FDR']:.3g}")
