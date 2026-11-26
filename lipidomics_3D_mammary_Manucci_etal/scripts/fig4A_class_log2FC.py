"""log2 fold change per lipid class (rBM vs nT), one point per biological replicate.

Data: data/lipidomics_cells_DNA_normalized.csv (unrounded, DNA-normalized).
Class sums: all species of each class; the duplicated feature pPE (p18:1/14:0)* is counted once.
Statistics: Student's t-test (equal var.) on log10 class sums, BH FDR over 20 classes.

Usage: python fig4A_class_log2FC.py <NormDNA.csv> <out_prefix> pontos|replicas
  pontos   -> original style (nT grey, rBM black)
  replicas -> points coloured by replicate number
"""
import sys, csv
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

CSV, OUT, MODE = sys.argv[1], sys.argv[2], sys.argv[3]
PPE_DUP = "pPE (p18:1/14:0)*"
ORDER = ["1G-cer", "Cer", "SM", "FFA", "ACar", "oPE", "PE", "pPE", "oPC", "pPC",
         "PC", "PG", "PI", "CL", "CoQ", "FC", "CE", "ADG", "DAG", "TAG"]
LABEL = {"1G-cer": "1G-Cer"}
CATEGORY = {"1G-cer": "Sphingolipids", "Cer": "Sphingolipids", "SM": "Sphingolipids",
            "FFA": "Fatty Acyls", "ACar": "Fatty Acyls", "CoQ": "Prenols", "FC": "Sterols", "CE": "Sterols",
            "ADG": "Glycerolipids", "DAG": "Glycerolipids", "TAG": "Glycerolipids"}
CATS = ["Sphingolipids", "Fatty Acyls", "Glycerophospholipids", "Prenols", "Sterols", "Glycerolipids"]
cat_col = dict(zip(CATS, plt.cm.viridis(np.linspace(0.0, 0.88, len(CATS)))))

def species_class(sp):
    c = sp.split(" ")[0]
    return "CoQ" if c.startswith("Q-") else c

rows = list(csv.reader(open(CSV, encoding="utf-8-sig")))
hdr = rows[0]; samples = {r[0]: r for r in rows[1:]}
nT = [f"nT{i}" for i in range(1, 5)]; rBM = [f"lrECM{i}" for i in range(1, 5)]

def class_sum(sample, c):
    r = samples[sample]
    return sum(float(v) for sp, v in zip(hdr[2:], r[2:]) if species_class(sp) == c and sp != PPE_DUP)

res = {}
for c in ORDER:
    a = np.array([class_sum(s, c) for s in nT]); b = np.array([class_sum(s, c) for s in rBM])
    ref = np.exp(np.mean(np.log(a)))                        # geometric mean of nT
    res[c] = dict(nT=np.log2(a / ref), rBM=np.log2(b / ref),
                  p=stats.ttest_ind(np.log10(b), np.log10(a)).pvalue)
ps = np.array([res[c]["p"] for c in ORDER]); m = len(ps); idx = np.argsort(ps)
q = np.minimum.accumulate((ps[idx] * m / np.arange(1, m + 1))[::-1])[::-1]
fdr = np.empty(m); fdr[idx] = np.minimum(q, 1)
for c, f in zip(ORDER, fdr): res[c]["fdr"] = f

mark = lambda f: "***" if f < .001 else "**" if f < .01 else "*" if f < .05 else "#" if f < .10 else ""
sem = lambda x: x.std(ddof=1) / np.sqrt(len(x))
REP_COL = ["#E69F00", "#D81B60", "#1E88E5", "#5D4037"]
REP_DY = [-0.09, -0.03, 0.03, 0.09]
Y_NT, Y_RBM = -0.22, 0.06

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig, ax = plt.subplots(figsize=(10.3, 11.4))
for y, c in enumerate(ORDER):
    d = res[c]; col = cat_col[CATEGORY.get(c, "Glycerophospholipids")]
    mb = d["rBM"].mean()
    ax.barh(y + Y_RBM, mb, height=0.6, color=col, alpha=0.35 if d["fdr"] >= 0.10 else 1.0, zorder=1, lw=0)
    ax.errorbar(mb, y + Y_RBM, xerr=sem(d["rBM"]), fmt="none", ecolor="black", elinewidth=1.2, capsize=3, zorder=3)
    ax.errorbar(d["nT"].mean(), y + Y_NT, xerr=sem(d["nT"]), fmt="none", ecolor="0.5", elinewidth=1.2,
                capsize=3, zorder=3)
    for k in range(4):
        if MODE == "replicas":
            ax.scatter(d["nT"][k], y + Y_NT + REP_DY[k], s=42, marker="o", facecolor="white",
                       edgecolor=REP_COL[k], linewidth=1.6, zorder=4)
            ax.scatter(d["rBM"][k], y + Y_RBM + REP_DY[k], s=48, marker="D", facecolor=REP_COL[k],
                       edgecolor="black", linewidth=0.6, zorder=5)
        else:
            ax.scatter(d["nT"][k], y + Y_NT + REP_DY[k], s=36, facecolor="0.68", edgecolor="0.35",
                       linewidth=0.7, zorder=4)
            ax.scatter(d["rBM"][k], y + Y_RBM + REP_DY[k], s=36, facecolor="0.2", edgecolor="0.2", zorder=5)
    s = mark(d["fdr"])
    if s:
        xmin = min(d["rBM"].min(), mb - sem(d["rBM"]), d["nT"].min())
        ax.text(xmin - 0.05, y + Y_RBM, s, ha="right", va="center", fontsize=12)
    if y: ax.axhline(y - 0.5, color="0.92", lw=0.8, zorder=0)

ax.axvline(0, color="black", lw=1, zorder=2)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([LABEL.get(c, c) for c in ORDER], fontsize=12)
ax.tick_params(axis="y", length=0)
ax.set_ylim(len(ORDER) - 0.5, -0.5); ax.set_xlim(-2.28, 1.1)
ax.set_xlabel(r"log$_2$ fold change (relative to nT)", fontsize=12)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)

leg_kw = dict(loc="upper left", frameon=False, alignment="left", fontsize=11, title_fontsize=11)
h3 = [Patch(color=cat_col[c], label=c) for c in CATS]
if MODE == "replicas":
    h1 = [Line2D([], [], ls="", marker="o", mfc="white", mec="0.3", mew=1.6, ms=7, label="nT"),
          Line2D([], [], ls="", marker="D", mfc="0.3", mec="black", mew=0.6, ms=7, label="rBM")]
    h2 = [Line2D([], [], ls="", marker="s", mfc=REP_COL[k], mec=REP_COL[k], ms=9, label=f"Replicate {k+1}")
          for k in range(4)]
    fig.legend(handles=h1, title="Condition (marker)", bbox_to_anchor=(0.765, 0.985), **leg_kw)
    fig.legend(handles=h2, title="Biological replicate (color)", bbox_to_anchor=(0.765, 0.905), **leg_kw)
    fig.legend(handles=h3, title="rBM bar: category", bbox_to_anchor=(0.765, 0.775), **leg_kw)
    first = "Each point is one biological replicate (n = 4 per condition; color = replicate number), expressed as\n"
else:
    h1 = [Line2D([], [], ls="", marker="o", mfc="0.68", mec="0.35", ms=7, label="nT"),
          Line2D([], [], ls="", marker="o", mfc="0.2", mec="0.2", ms=7, label="rBM")]
    fig.legend(handles=h1, title="Samples", bbox_to_anchor=(0.765, 0.985), **leg_kw)
    fig.legend(handles=h3, title="rBM bar: category", bbox_to_anchor=(0.765, 0.905), **leg_kw)
    first = "Each point is one sample (n = 4 per condition), expressed as\n"

fig.subplots_adjust(left=0.08, right=0.75, top=0.985, bottom=0.10)
fig.text(0.015, 0.012, first +
         "log$_2$(value / geometric mean of nT). Bars (rBM) and error bars: mean ± SEM. Faded bars: FDR ≥ 0.10.\n"
         "Student's t-test on log$_{10}$ data, Benjamini–Hochberg FDR across 20 classes: "
         "*** < 0.001; ** < 0.01; * < 0.05; # < 0.10.", fontsize=9.5, color="0.3", va="bottom")
fig.savefig(OUT + ".pdf")
if MODE == "replicas": fig.savefig(OUT + ".png", dpi=200)
print(MODE, {c: mark(res[c]["fdr"]) for c in ORDER if mark(res[c]["fdr"])})
