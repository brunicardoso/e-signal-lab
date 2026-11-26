"""Fig. 4B - acyl-chain remodelling (left/A: composition, % of acyl chains in class; right/B: abundance normalized to DNA).

Data: data/lipidomics_cells_DNA_normalized.csv (unrounded, DNA-normalized) and the
species -> acyl-chain map of data/fatty_acid_species_map.xlsx (sheet 'Mapa espécies').
Calculations identical to Supplementary Table S4:
  abundance of FA x in class C = sum(species value x number of x chains); % = share within class;
  log2 FC = log2(mean rBM / mean nT); t-test on log10 abundance (B) or on % (A); BH FDR per panel.

Usage: python fig4B_fa_remodeling.py <NormDNA.csv> <fatty_acid_species_map.xlsx> <out.pdf>
"""
import sys, csv, re
import numpy as np
from scipy import stats
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

CSV, XLSX, OUT = sys.argv[1:4]
NT = [f"nT{i}" for i in range(1, 5)]; RB = [f"lrECM{i}" for i in range(1, 5)]

rows = list(csv.reader(open(CSV, encoding="utf-8-sig")))
hdr = rows[0]
data = {r[0]: dict(zip(hdr[2:], map(float, r[2:]))) for r in rows[1:]}
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
mp = list(wb["Mapa espécies"].iter_rows(values_only=True))
FAS = list(mp[0][7:])
smap = {r[0]: (r[1], {fa: (r[7 + i] or 0) for i, fa in enumerate(FAS)}) for r in mp[1:]}

COLS = [("Glycerophospholipids", ["PC", "oPC", "pPC", "PE", "oPE", "pPE", "PG", "PI", "CL"]),
        ("Sphingolipids", ["Cer", "1G-cer", "SM"]),
        ("Glycerolipids", ["DAG", "ADG", "TAG"]),
        ("Ster.", ["CE"]),
        ("FA-yl", ["FFA", "ACar"])]
COLLAB = {"CL": "CL†", "1G-cer": "1G-Cer"}
ROWS = [("SFA", ["12:0", "14:0", "15:0", "16:0", "17:0", "18:0", "19:0", "20:0", "22:0", "24:0"]),
        ("MUFA", ["16:1", "17:1", "18:1", "20:1", "22:1", "24:1"]),
        ("PUFA", ["18:2", "18:3", "20:2", "20:3", "20:4", "20:5", "22:2", "22:3", "22:4", "22:5", "22:6", "24:2"])]
ROWLAB = {"20:4": "20:4 (ARA)", "20:5": "20:5 (EPA)", "22:5": "22:5 (DPA)", "22:6": "22:6 (DHA)"}

classes = [c for _, cs in COLS for c in cs]
# abundance per (class, FA) per sample
ab = {}
for c in classes:
    for fa in FAS:
        v = {s: sum(cnt[fa] * data[s][sp] for sp, (cl, cnt) in smap.items() if cl == c) for s in NT + RB}
        if any(v.values()): ab[(c, fa)] = v
pct = {k: {s: 100 * ab[k][s] / sum(ab[j][s] for j in ab if j[0] == k[0]) for s in NT + RB} for k in ab}

def bh(p):
    p = np.asarray(p); m = len(p); i = np.argsort(p); q = p[i] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; o = np.empty(m); o[i] = np.minimum(q, 1); return o

def panel_stats(table, log_test):
    keys = list(table); l2, p = {}, []
    for k in keys:
        a = np.array([table[k][s] for s in NT]); b = np.array([table[k][s] for s in RB])
        l2[k] = np.log2(b.mean() / a.mean())
        p.append(stats.ttest_ind(np.log10(b), np.log10(a)).pvalue if log_test else stats.ttest_ind(b, a).pvalue)
    return l2, dict(zip(keys, bh(p)))

# row order of the original workbook defines the FDR families: all (class, FA) rows present in the data
mark = lambda q: "**" if q < 0.01 else "*" if q < 0.05 else "#" if q < 0.10 else ""
minor = {k: np.mean([pct[k][s] for s in NT]) < 1 for k in pct}

def draw(ax, cax, l2, q, cap, title, letter):
    cmap = plt.cm.RdBu_r; norm = TwoSlopeNorm(0, -cap, cap)
    x = 0; xpos = []
    for gi, (gname, cs) in enumerate(COLS):
        for c in cs: xpos.append((c, x)); x += 1
        x += 0.35
    y = 0; ypos = []
    for ri, (rname, fas) in enumerate(ROWS):
        for fa in fas: ypos.append((fa, y)); y += 1
        y += 0.35
    for c, xx in xpos:
        for fa, yy in ypos:
            k = (c, fa)
            if k not in l2:
                ax.add_patch(plt.Rectangle((xx, yy), 0.94, 0.94, color="#f2f2f2", lw=0)); continue
            v = l2[k]
            ax.add_patch(plt.Rectangle((xx, yy), 0.94, 0.94, color=cmap(norm(np.clip(v, -cap, cap))), lw=0))
            m = mark(q[k])
            if m:
                ax.text(xx + 0.47, yy + 0.5, m, ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > 0.7 * cap else "black")
            if minor[k]:
                ax.plot(xx + 0.86, yy + 0.12, marker="o", ms=1.6, color="black")
    ax.set_xlim(-0.1, x); ax.set_ylim(y, -0.1)
    ax.set_xticks([xx + 0.47 for _, xx in xpos]); ax.set_xticklabels([COLLAB.get(c, c) for c, _ in xpos], rotation=90)
    ax.xaxis.tick_top(); ax.tick_params(length=0, labelsize=8)
    ax.set_yticks([yy + 0.47 for _, yy in ypos]); ax.set_yticklabels([ROWLAB.get(f, f) for f, _ in ypos], fontsize=8)
    for s in ax.spines.values(): s.set_visible(False)
    # group brackets
    for gname, cs in COLS:
        xs = [xx for c, xx in xpos if c in cs]
        ax.plot([xs[0], xs[-1] + 0.94], [y + 0.2] * 2, color="0.5", lw=0.8, clip_on=False)
        ax.text((xs[0] + xs[-1] + 0.94) / 2, y + 0.75, gname, ha="center", va="top", fontsize=7, color="0.3")
    for rname, fas in ROWS:
        ys = [yy for f, yy in ypos if f in fas]
        ax.plot([-2.6, -2.6], [ys[0], ys[-1] + 0.94], color="0.5", lw=0.8, clip_on=False)
        ax.text(-3.0, (ys[0] + ys[-1] + 0.94) / 2, rname, rotation=90, ha="right", va="center", fontsize=8, color="0.3")
    ax.set_title(title, fontsize=9, pad=48)
    ax.text(-0.12, 1.12, letter, transform=ax.transAxes, fontsize=14, fontweight="bold")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = plt.colorbar(sm, cax=cax); cb.set_label("log$_2$ FC (rBM / nT)", fontsize=7); cb.ax.tick_params(labelsize=7)

fig = plt.figure(figsize=(14.5, 9.3))
axA = fig.add_axes([0.065, 0.17, 0.37, 0.70]); caxA = fig.add_axes([0.445, 0.45, 0.008, 0.25])
axB = fig.add_axes([0.565, 0.17, 0.37, 0.70]); caxB = fig.add_axes([0.945, 0.45, 0.008, 0.25])
l2A, qA = panel_stats(pct, False); l2B, qB = panel_stats(ab, True)
draw(axA, caxA, l2A, qA, 1.0, "Composition (% of acyl chains in class)", "A")
draw(axB, caxB, l2B, qB, 1.5, "Abundance (normalized to DNA)", "B")
fig.text(0.02, 0.015,
         "Each cell: log$_2$ ratio (rBM/nT) of a fatty acid esterified in a lipid class (n = 4). (A) share of the fatty acid in the acyl chains of the class; Student's t-test on percentages; colour scale capped at ±1.\n"
         "(B) abundance of the fatty acid in the class (Σ species × number of chains, normalized to DNA); Student's t-test on log$_{10}$ abundance; capped at ±1.5. Benjamini–Hochberg FDR across all\n"
         "class × fatty-acid combinations in each panel: ** FDR < 0.01; * < 0.05; # < 0.10. Light grey: not detected in the class. Small dot: minor component (< 1% of the class in nT).\n"
         "Ether (alkyl/alkenyl) chains and sphingoid bases are not included. † CL: high technical variability in QC samples (CV > 40%); interpret with caution. Ster., sterol esters; FA-yl, fatty acyls.",
         fontsize=7, color="0.3", va="bottom")
fig.savefig(OUT)
nA = sum(bool(mark(v)) for v in qA.values()); nB = sum(bool(mark(v)) for v in qB.values())
print(f"saved {OUT}; marked cells A={nA} B={nB}; tests A={len(qA)} B={len(qB)}")
