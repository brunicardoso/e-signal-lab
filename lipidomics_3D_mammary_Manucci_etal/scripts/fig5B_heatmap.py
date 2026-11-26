"""Fig. 5B - clustered heatmap of the significant species .

Significant species: FDR < 0.10 and |log2FC| > log2(1.5), read from the volcano statistics CSV
(fig5A_volcano.py output). Rows are z-scored DNA-normalized abundances; seaborn clustermap
defaults (Euclidean distance, average linkage) for rows and columns.
Display labels: '1G-cer' is shown as '1G-Cer' (data names unchanged); 'lrECM' samples are shown as 'rBM'.

Usage: python fig5B_heatmap.py <NormDNA.csv> <volcano_stats.csv> <out.png> [--raw-labels]
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

CSV, STATS, OUT = sys.argv[1:4]
RAW = "--raw-labels" in sys.argv
GROUP_NUM, GROUP_DEN = "lrECM", "nT"
COLOR_rBM, COLOR_nT = "#8BC34A", "#9E9E9E"
label = (lambda s: s) if RAW else (lambda s: s.replace("1G-cer", "1G-Cer"))

df = pd.read_csv(CSV)
st = pd.read_csv(STATS)
sig = set(st.loc[st["Significant"].astype(str) == "True", "Identifier"])
species = [c for c in df.columns if c not in ("Sample", "Group")]
sig_species = [c for c in species if c in sig]          # same order as the notebook (column order)

d = df[df["Group"].isin([GROUP_NUM, GROUP_DEN])].copy()
d["Sample"] = d["Sample"].str.replace("lrECM", "rBM", regex=False)
col_colors = d.set_index("Sample")["Group"].map({GROUP_NUM: COLOR_rBM, GROUP_DEN: COLOR_nT})
heat = d.set_index("Sample")[sig_species].T.apply(pd.to_numeric, errors="coerce")
heat.index = [label(s) for s in heat.index]               # display only; values and order unchanged

g = sns.clustermap(heat, z_score=0, cmap="vlag", figsize=(12, max(8, len(sig_species) * 0.3)),
                   col_colors=col_colors, col_cluster=True, row_cluster=True, linewidths=0.5,
                   cbar_pos=(0.02, 0.8, 0.03, 0.18))
g.ax_cbar.set_ylabel("Z-Score", fontsize=12, fontweight="bold")
g.ax_heatmap.set_xlabel("Samples", fontsize=14, fontweight="bold")
g.ax_heatmap.set_ylabel("Significant Lipids", fontsize=14, fontweight="bold")
g.ax_col_dendrogram.legend(handles=[Patch(facecolor=COLOR_nT, label="nT"), Patch(facecolor=COLOR_rBM, label="rBM")],
                           loc="upper right", title="Groups")
plt.savefig(OUT, dpi=300, bbox_inches="tight")
print("saved", OUT, "| rows:", [t.get_text() for t in g.ax_heatmap.get_yticklabels()],
      "| cols:", [t.get_text() for t in g.ax_heatmap.get_xticklabels()])
