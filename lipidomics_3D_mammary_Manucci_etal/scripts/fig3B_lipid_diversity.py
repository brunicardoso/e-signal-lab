"""Fig. 3B - number of lipid species identified per class in the 3D structures, coloured by category.

Data: data/lipidomics_cells_DNA_normalized.csv (species = columns). Every annotated species is counted,
including both names of the ambiguous feature pPE (p18:1/14:0)* / pPE (p14/18:1)* (211 species).
Colours: same category palette as Fig. 4A (fig4A_class_log2FC.py).

Usage: python fig3B_lipid_diversity.py <NormDNA.csv> <out_prefix>   -> <out_prefix>.png/.pdf/.svg
"""
import sys, csv
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CSV, OUT = sys.argv[1], sys.argv[2]
ORDER = ["1G-cer", "Cer", "SM", "FFA", "ACar", "CL", "oPC", "pPC", "PC", "oPE", "pPE", "PE", "PG", "PI",
         "CoQ", "CE", "FC", "DAG", "ADG", "TAG"]
LABEL = {"1G-cer": "1G-Cer"}
CATEGORY = {"1G-cer": "Sphingolipids", "Cer": "Sphingolipids", "SM": "Sphingolipids",
            "FFA": "Fatty Acyls", "ACar": "Fatty Acyls", "CoQ": "Prenols", "FC": "Sterols", "CE": "Sterols",
            "ADG": "Glycerolipids", "DAG": "Glycerolipids", "TAG": "Glycerolipids"}
CATS = ["Sphingolipids", "Fatty Acyls", "Glycerophospholipids", "Prenols", "Sterols", "Glycerolipids"]
cat_col = dict(zip(CATS, plt.cm.viridis(np.linspace(0.0, 0.88, len(CATS)))))

def species_class(sp):
    c = sp.split(" ")[0]
    return "CoQ" if c.startswith("Q-") else c

hdr = next(csv.reader(open(CSV, encoding="utf-8-sig")))
n = Counter(species_class(sp) for sp in hdr[2:])
assert set(n) == set(ORDER), set(n) ^ set(ORDER)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 15})
fig, ax = plt.subplots(figsize=(11.6, 4.6))
x = np.arange(len(ORDER))
cols = [cat_col[CATEGORY.get(c, "Glycerophospholipids")] for c in ORDER]
ax.bar(x, [n[c] for c in ORDER], width=0.78, color=cols, edgecolor="white", linewidth=0.8, zorder=2)
ax.set_xticks(x); ax.set_xticklabels([LABEL.get(c, c) for c in ORDER], rotation=45, ha="right",
                                     rotation_mode="anchor", fontsize=15)
ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.6, len(ORDER) - 0.4); ax.set_ylim(0, 60); ax.set_yticks(range(0, 61, 20))
ax.tick_params(axis="y", labelsize=15)
ax.set_ylabel("Total lipid diversity\n(number of identified species)", fontsize=15)
ax.yaxis.grid(True, color="0.9", lw=0.8, zorder=0)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(handles=[Patch(color=cat_col[c], label=c) for c in CATS], loc="upper left",
          bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=15)
fig.subplots_adjust(left=0.115, right=0.70, top=0.95, bottom=0.22)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.1)
print({c: n[c] for c in ORDER}, "total", sum(n.values()),
      {k: sum(n[c] for c in ORDER if CATEGORY.get(c, "Glycerophospholipids") == k) for k in CATS})
