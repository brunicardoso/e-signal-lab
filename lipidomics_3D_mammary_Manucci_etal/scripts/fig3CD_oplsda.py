"""Fig. 3C-D - OPLS-DA of the cell lipidome (nT vs rBM): score plot (3C) and S-plot (3D).

Python re-implementation of the MetaboAnalyst analysis that produced the original figure
(exported as data/metaboanalyst_oplsda_loadings.csv). The settings below were recovered by matching
that export and the published R2X values (28.8 % / 24.5 %):

  * input: DNA-normalized concentrations (data/lipidomics_cells_DNA_normalized.csv); pooled QC
    samples (Group == "Pool") excluded -> 8 samples (nT n=4, rBM n=4), 211 lipids
  * log2 transform, then auto-scaling (mean-centred, divided by SD; computed on the 8 samples)
  * OPLS-DA (NIPALS; Trygg & Wold 2002, as in ropls/MetaboAnalyst): 1 predictive +
    1 orthogonal component; y = class dummy, mean-centred
  * S-plot (as in MetaboAnalyst PlotOPLS.Splot):
        p[1]       = cov(log2 x_j, t_p)                      (x axis, "covariance")
        p(corr)[1] = cor(x_j, t_p)                           (y axis, "correlation")
    where t_p is the predictive score. MetaboAnalyst additionally adds jitter() to p[1];
    that random jitter is not reproduced here.
  * sign of the predictive component fixed so that nT is positive (as in the published figure)

Model validation (not given in MetaboAnalyst's plots, added here):
  * R2X (predictive, orthogonal), R2Y, Q2 by 7-fold cross-validation (ropls segment scheme:
    sample i goes to fold i mod 7, so with n = 8 one fold holds 2 samples)
  * exact permutation test: all C(8,4) = 70 relabellings of the 8 samples

Usage: python fig3CD_oplsda.py <NormDNA.csv> <out_prefix> [MetaboAnalyst_export.csv]
  -> <out_prefix>_3C_scores.{png,pdf,svg}, <out_prefix>_3D_splot.{png,pdf,svg},
     <out_prefix>_scores.csv, <out_prefix>_splot.csv, <out_prefix>_model.txt
"""
import sys, csv, itertools, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CSV, OUT = sys.argv[1], sys.argv[2]
REF = sys.argv[3] if len(sys.argv) > 3 else None
GROUP_POS, GROUP_NEG, NEG_LABEL = "nT", "lrECM", "rBM"
N_ORTHO, CV_FOLDS = 1, 7

# same category palette as Fig. 3B (plot_lipid_diversity.py)
CATEGORY = {"1G-cer": "Sphingolipids", "Cer": "Sphingolipids", "SM": "Sphingolipids",
            "FFA": "Fatty Acyls", "ACar": "Fatty Acyls", "CoQ": "Prenols", "FC": "Sterols", "CE": "Sterols",
            "ADG": "Glycerolipids", "DAG": "Glycerolipids", "TAG": "Glycerolipids"}
CATS = ["Sphingolipids", "Fatty Acyls", "Glycerophospholipids", "Prenols", "Sterols", "Glycerolipids"]
cat_col = dict(zip(CATS, plt.cm.viridis(np.linspace(0.0, 0.88, len(CATS)))))
GROUP_COL = {"nT": plt.cm.viridis(0.33), "rBM": plt.cm.viridis(0.83)}

def species_class(sp):
    c = sp.split(" ")[0]
    return "CoQ" if c.startswith("Q-") else c

# ---------------------------------------------------------------- data
rows = list(csv.reader(open(CSV, encoding="utf-8-sig")))
hdr, body = rows[0], [r for r in rows[1:] if r[1] in (GROUP_POS, GROUP_NEG)]
species = hdr[2:]
samples = [r[0] for r in body]
groups = [GROUP_POS if r[1] == GROUP_POS else NEG_LABEL for r in body]
raw = np.array([[float(v) for v in r[2:]] for r in body])
assert (raw > 0).all()
L2 = np.log2(raw)
X = (L2 - L2.mean(0)) / L2.std(0, ddof=1)               # auto-scaling
y = np.array([1.0 if g == GROUP_POS else 0.0 for g in groups])
n, m = X.shape

# ---------------------------------------------------------------- OPLS (single y, NIPALS)
def opls_fit(X, y, n_ortho):
    """X, y already centred. Returns predictive/orthogonal weights, loadings, scores."""
    Xr = X.copy()
    Wo, Po, To = [], [], []
    for _ in range(n_ortho):
        w = Xr.T @ y / (y @ y); w /= np.linalg.norm(w)
        t = Xr @ w
        p = Xr.T @ t / (t @ t)
        wo = p - (w @ p) * w; wo /= np.linalg.norm(wo)
        to = Xr @ wo
        po = Xr.T @ to / (to @ to)
        Xr = Xr - np.outer(to, po)
        Wo.append(wo); Po.append(po); To.append(to)
    w = Xr.T @ y / (y @ y); w /= np.linalg.norm(w)
    t = Xr @ w
    p = Xr.T @ t / (t @ t)
    c = (y @ t) / (t @ t)
    return dict(w=w, p=p, t=t, c=c, Wo=np.array(Wo), Po=np.array(Po), To=np.array(To).T)

def opls_predict(model, Xnew):
    Xr = Xnew.copy()
    for wo, po in zip(model["Wo"], model["Po"]):
        to = Xr @ wo
        Xr = Xr - np.outer(to, po)
    return (Xr @ model["w"]) * model["c"]

def r2y_q2(X, y, n_ortho, folds=CV_FOLDS):
    yc = y - y.mean()
    mod = opls_fit(X - X.mean(0), yc, n_ortho)
    r2y = 1 - ((yc - mod["t"] * mod["c"]) ** 2).sum() / (yc ** 2).sum()
    seg = np.arange(len(y)) % min(folds, len(y))
    press = 0.0
    for k in np.unique(seg):
        tr, te = seg != k, seg == k
        mx, my = X[tr].mean(0), y[tr].mean()
        mk = opls_fit(X[tr] - mx, y[tr] - my, n_ortho)
        press += ((y[te] - (opls_predict(mk, X[te] - mx) + my)) ** 2).sum()
    q2 = 1 - press / (yc ** 2).sum()
    return r2y, q2

Xc = X - X.mean(0)
ss_x = (Xc ** 2).sum()
mod = opls_fit(Xc, y - y.mean(), N_ORTHO)
t_p, t_o = mod["t"], mod["To"][:, 0]
sign = 1.0 if t_p[y == 1].mean() > 0 else -1.0            # nT positive
t_p, p_p = sign * t_p, sign * mod["p"]
r2x_p = (np.outer(t_p, p_p) ** 2).sum() / ss_x
r2x_o = (np.outer(mod["To"][:, 0], mod["Po"][0]) ** 2).sum() / ss_x
r2y, q2 = r2y_q2(X, y, N_ORTHO)

# exact permutation test over all 70 relabellings (includes the observed labelling)
perm = []
for pos in itertools.combinations(range(n), int(y.sum())):
    yp = np.zeros(n); yp[list(pos)] = 1
    perm.append(r2y_q2(X, yp, N_ORTHO))
perm = np.array(perm)
p_r2y = (perm[:, 0] >= r2y - 1e-12).mean()
p_q2 = (perm[:, 1] >= q2 - 1e-12).mean()

# Q2 for other numbers of orthogonal components (for reference)
q2_by_north = {k: r2y_q2(X, y, k) for k in range(0, 4)}

# S-plot
L2c = L2 - L2.mean(0)
p1 = L2c.T @ t_p / (n - 1)                                     # cov(log2 x_j, t_p)
pcorr = p1 / (L2.std(0, ddof=1) * t_p.std(ddof=1))            # cor(x_j, t_p)

# ---------------------------------------------------------------- outputs (tables)
with open(f"{OUT}_scores.csv", "w", newline="") as f:
    w_ = csv.writer(f); w_.writerow(["Sample", "Group", "t_pred", "t_ortho"])
    for s, g, a, b in zip(samples, groups, t_p, t_o): w_.writerow([s, g, f"{a:.4f}", f"{b:.4f}"])
with open(f"{OUT}_splot.csv", "w", newline="") as f:
    w_ = csv.writer(f); w_.writerow(["Lipid", "Class", "Category", "p[1]", "p(corr)[1]"])
    for sp, a, b in zip(species, p1, pcorr):
        cl = species_class(sp)
        w_.writerow([sp, cl, CATEGORY.get(cl, "Glycerophospholipids"), f"{a:.5f}", f"{b:.5f}"])

lines = [f"OPLS-DA nT vs rBM: {n} samples x {m} lipids; log2 + auto-scaling; 1 predictive + {N_ORTHO} orthogonal",
         f"R2X predictive  = {100*r2x_p:.1f} %",
         f"R2X orthogonal  = {100*r2x_o:.1f} %",
         f"R2Y             = {r2y:.3f}",
         f"Q2 ({CV_FOLDS}-fold CV)  = {q2:.3f}",
         f"Permutation test (all {len(perm)} relabellings): p(R2Y) = {p_r2y:.3f}, p(Q2) = {p_q2:.3f}",
         "Q2 by number of orthogonal components: " +
         ", ".join(f"{k}: R2Y={a:.3f} Q2={b:.3f}" for k, (a, b) in q2_by_north.items())]

# optional check against the MetaboAnalyst export
if REF:
    norm = lambda s: re.sub(r"[^A-Za-z0-9]", "", s).lower()
    ref = {norm(r[0]): (float(r[1]), float(r[2])) for r in list(csv.reader(open(REF, encoding="utf-8-sig")))[1:]}
    alias = {norm("oPE (o16/20:1)"): norm("pPE o16/201")}   # species renamed after the MetaboAnalyst run
    d1, d2 = [], []
    for sp, a, b in zip(species, p1, pcorr):
        k = norm(sp); k = k if k in ref else alias.get(k, k)
        d1.append(abs(a - ref[k][0])); d2.append(abs(b - ref[k][1]))
    d1, d2 = np.array(d1), np.array(d2)
    lines.append(f"vs MetaboAnalyst export: |d p[1]| median {np.median(d1):.5f}, max {d1.max():.4f}; "
                 f"|d p(corr)| median {np.median(d2):.6f}, max {d2.max():.4f}; "
                 f"lipids with |d p(corr)| > 1e-4: {(d2 > 1e-4).sum()}")
open(f"{OUT}_model.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------------------------------------------------------------- Fig. 3C - score plot
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 15})
fig, ax = plt.subplots(figsize=(5.6, 4.6))
for g in ("nT", "rBM"):
    k = np.array(groups) == g
    ax.scatter(t_p[k], t_o[k], s=110, color=GROUP_COL[g], edgecolor="white", linewidth=1.2, zorder=3, label=g)
lim_x = np.ceil(np.abs(t_p).max() / 5) * 5; lim_y = np.ceil(np.abs(t_o).max() / 10) * 10
ax.set_xlim(-lim_x, lim_x); ax.set_ylim(-lim_y, lim_y)
ax.set_yticks(np.arange(-lim_y, lim_y + 1, 5))
ax.axhline(0, color="0.85", lw=0.8, zorder=0); ax.axvline(0, color="0.85", lw=0.8, zorder=0)
ax.set_xlabel(f"T score [1] ({100*r2x_p:.1f}%)")
ax.set_ylabel(f"Orthogonal T score [1] ({100*r2x_o:.1f}%)")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, handletextpad=0.2)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}_3C_scores.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.1)

# ---------------------------------------------------------------- Fig. 3D - S-plot
# cross-shaped axes through the origin (as in the original MetaboAnalyst/Prism panel)
FS = 19.5                                   # 1.3 x the 15 pt used in Fig. 3B
AX_W, AX_H = 4.56, 3.54                     # plot area in inches (30 % narrower than the first version)
fig = plt.figure(figsize=(AX_W + 7.0, AX_H + 1.6))
W, H = fig.get_size_inches()
ax = fig.add_axes([1.3 / W, 1.1 / H, AX_W / W, AX_H / H])
cats = np.array([CATEGORY.get(species_class(sp), "Glycerophospholipids") for sp in species])
for c in CATS:
    k = cats == c
    ax.scatter(p1[k], pcorr[k], s=45, color=cat_col[c], alpha=0.8, edgecolor="white", linewidth=0.5, zorder=3)
lim = np.ceil(np.abs(p1).max() / 2.5) * 2.5
ax.set_xlim(-lim, lim); ax.set_ylim(-1.05, 1.05)
for sp in ("left", "bottom"):
    ax.spines[sp].set_position("zero"); ax.spines[sp].set_linewidth(1.5); ax.spines[sp].set_zorder(2)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
xt = [t for t in np.arange(-np.floor(lim / 5) * 5, lim, 5) if t != 0]
yt = [-1, -0.5, 0.5, 1]
ax.set_xticks(xt); ax.set_yticks(yt)
ax.tick_params(direction="inout", width=1.5, length=9, labelbottom=False, labelleft=False)
# tick labels placed in the empty quadrants (upper-left / lower-right) so they never cover points
dx, dy = 0.025 * (2 * lim), 0.05
for t in xt:
    ax.text(t, dy if t < 0 else -dy, f"{t:g}".replace("-", "\u2212"), ha="center", va="bottom" if t < 0 else "top", fontsize=FS)
for t in yt:
    ax.text(-dx if t > 0 else dx, t, f"{t:g}".replace("-", "\u2212"), ha="right" if t > 0 else "left", va="center", fontsize=FS)
ax.set_xlabel("p[1] (covariance)", fontsize=FS)
ax.set_ylabel("p(corr)[1]\n(correlation)", fontsize=FS, linespacing=1.0)
ax.xaxis.set_label_coords(0.5, -0.06, transform=ax.transAxes); ax.xaxis.label.set_va("top")
ax.yaxis.set_label_coords(-0.06, 0.5, transform=ax.transAxes); ax.yaxis.label.set_va("bottom")
ax.legend(handles=[Line2D([], [], ls="", marker="o", ms=11, color=cat_col[c], label=c) for c in CATS],
          loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, handletextpad=0.2,
          fontsize=FS, labelspacing=0.35)
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"{OUT}_3D_splot.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.1)
