# -*- coding: utf-8 -*-
"""Step 5: figures + analysis summary for the cholesterol-metabolism signature study."""
import os, json
os.environ["LOKY_MAX_CPU_COUNT"] = "4"
from pathlib import Path as _P
(_P(r".\work\mplcache_chol").mkdir(parents=True, exist_ok=True))
os.environ["MPLCONFIGDIR"] = r".\work\mplcache_chol"
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.metrics import roc_auc_score

OUT = Path(r".\outputs\chol_metab_signature")
FIG = OUT / "figures"
FIG.mkdir(exist_ok=True)
RNG = 20260804

X = np.load(OUT / "tcga_Xlog.npy")
genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
tr = pd.read_csv(OUT / "tcga_traits.csv")
mb = pd.read_csv(OUT / "mb_traits.csv")
idx = {g: i for i, g in enumerate(genes)}

sf = pd.read_csv(OUT / "soft_threshold.csv")
mods = pd.read_csv(OUT / "modules_summary.csv")
mt = pd.read_csv(OUT / "module_trait.csv")
hub = pd.read_csv(OUT / "hub_genes.csv")
ml = json.load(open(OUT / "ml_summary.json", encoding="utf-8"))
pac = json.load(open(OUT / "consensus_summary.json", encoding="utf-8"))
sub = pd.read_csv(OUT / "tcga_subtypes.csv")
imm = pd.read_csv(OUT / "immune_scores.csv")
sub_stats = pd.read_csv(OUT / "subtype_stats.csv", index_col=0)

wgcna = json.load(open(OUT / "wgcna_summary.json", encoding="utf-8"))
sig_t = np.load(OUT / "signature_tcga.npy")
sig_m = np.load(OUT / "signature_mb.npy")

# ---------- Fig 1: soft threshold + module sizes ----------
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
ax[0].plot(sf.beta, sf.scale_free_R2, "o-", color="#2c7fb8")
ax[0].axhline(0.8, ls="--", color="grey", lw=0.8)
ax[0].axvline(wgcna["beta"], ls=":", color="red", lw=1)
ax[0].set_xlabel("Soft-threshold power (beta)"); ax[0].set_ylabel("Scale-free fit index R2")
ax[0].set_title("Soft-threshold selection (chosen beta = %d)" % wgcna["beta"])
mods2 = mods.sort_values("n_genes")
ax[1].bar(range(len(mods2)), mods2.n_genes, color="#74add1")
ax[1].set_xticks(range(len(mods2))); ax[1].set_xticklabels(mods2.module, rotation=90, fontsize=7)
ax[1].set_ylabel("Genes per module"); ax[1].set_title("%d WGCNA modules" % len(mods2))
fig.tight_layout(); fig.savefig(FIG / "fig1_wgcna.png", dpi=200); plt.close(fig)

# ---------- Fig 2: module-trait heatmap ----------
piv = mt.pivot(index="module", columns="trait", values="cor")
pv = mt.pivot(index="module", columns="trait", values="padj")
order = mt[mt.trait == "event"].sort_values("p").module.tolist()
order = order + [m for m in piv.index if m not in order]
piv = piv.loc[order, ["event", "time", "ER", "stage"]]
pv = pv.loc[order, ["event", "time", "ER", "stage"]]
fig, ax = plt.subplots(figsize=(5.4, 6.2))
im = ax.imshow(piv.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(piv.columns)
ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=7)
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        v = pv.values[i, j]
        ax.text(j, i, "%.2f" % piv.values[i, j], ha="center", va="center", fontsize=6,
                color="white" if abs(piv.values[i, j]) > 0.55 else "black")
        if v < 0.05:
            ax.text(j, i + 0.32, "*", ha="center", fontsize=7, color="red")
ax.set_title("Module-trait correlations (BH-FDR *)")
fig.colorbar(im, ax=ax, shrink=0.7)
fig.tight_layout(); fig.savefig(FIG / "fig2_module_trait.png", dpi=200); plt.close(fig)

# ---------- Fig 3: hub coefficients + KM ----------
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
hub2 = hub.sort_values("coef")
ax[0].barh(hub2.gene, hub2.coef, color=["#d62728" if c > 0 else "#2c7fb8" for c in hub2.coef])
ax[0].axvline(0, color="grey", lw=0.8)
ax[0].set_xlabel("LASSO coefficient"); ax[0].set_title("Hub gene coefficients (ER+ vs ER-)")

for a, (name, t, e, s, col) in enumerate([("TCGA-BRCA", tr.time, tr.event, sig_t, "#2c7fb8"),
                                           ("METABRIC", mb.time, mb.event, sig_m, "#d62728")], start=1):
    med = np.median(s)
    g1 = s <= med; g2 = s > med
    kmf1 = KaplanMeierFitter().fit(t[g1], e[g1], label="Low signature")
    kmf2 = KaplanMeierFitter().fit(t[g2], e[g2], label="High signature")
    lr = logrank_test(t[g1], t[g2], e[g1], e[g2])
    ax[a].plot(kmf1.survival_function_.index / 365.25, kmf1.survival_function_.values, color="#2c7fb8")
    ax[a].plot(kmf2.survival_function_.index / 365.25, kmf2.survival_function_.values, color="#d62728")
    ax[a].set_title("%s (log-rank p = %.3f)" % (name, lr.p_value))
    ax[a].set_xlabel("Years"); ax[a].set_ylabel("Survival probability")
fig.tight_layout(); fig.savefig(FIG / "fig3_signature.png", dpi=200); plt.close(fig)

# ---------- Fig 4: ROC + ER boxplot ----------
roc = pd.read_csv(OUT / "roc_metabric.csv")
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
ax[0].plot(roc.fpr, roc.tpr, color="#d62728", lw=2)
ax[0].plot([0, 1], [0, 1], ls="--", color="grey", lw=0.8)
ax[0].set_xlabel("1 - specificity"); ax[0].set_ylabel("Sensitivity")
ax[0].set_title("ER-status classifier in METABRIC (AUC = %.3f)" % ml["external_auc_metabric"])
for a, (name, s, e) in enumerate([("TCGA", sig_t, tr.ER.values), ("METABRIC", sig_m, mb.ER.values)], start=0):
    pos = s[e == 1]; neg = s[e == 0]
    bp = ax[a].boxplot([pos, neg], tick_labels=["ER+", "ER-"], patch_artist=True,
                       medianprops=dict(color="black"))
    bp["boxes"][0].set_facecolor("#d62728"); bp["boxes"][1].set_facecolor("#2c7fb8")
    ax[a].set_title(name); ax[a].set_ylabel("Signature score")
fig.tight_layout(); fig.savefig(FIG / "fig4_roc_er.png", dpi=200); plt.close(fig)

# ---------- Fig 5: consensus + subtype KM + immune + hub heatmap ----------
fig = plt.figure(figsize=(13.5, 9))
gs = fig.add_gridspec(2, 3)
ax = fig.add_subplot(gs[0, 0])
ks = sorted(pac["pac"])
ax.plot(ks, [pac["pac"][k] for k in ks], "o-", color="#2c7fb8")
ax.set_xlabel("k"); ax.set_ylabel("PAC"); ax.set_title("Consensus clustering (chosen k = %d)" % pac["k_best"])

ax = fig.add_subplot(gs[0, 1:])
lab = sub.subtype.values
for c in np.unique(lab):
    m = lab == c
    kmf = KaplanMeierFitter().fit(sub.time[m], sub.event[m], label="C%d (n=%d)" % (c, m.sum()))
    kmf.plot(ax=ax, ci_show=False)
ax.set_title("PFI by cholesterol-metabolism subtype (multivariate log-rank p = %.3f)" %
             json.load(open(OUT / "subtype_summary.json", encoding="utf-8"))["multivariate_logrank_p"])
ax.set_xlabel("Days"); ax.set_ylabel("PFI probability")

ax = fig.add_subplot(gs[1, 0])
imm2 = imm.set_index("subtype").sort_index()
colz = [c for c in imm2.columns if c != "subtype"]
im = ax.imshow(imm2[colz].T.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
ax.set_yticks(range(len(colz))); ax.set_yticklabels(colz, fontsize=8)
ax.set_xticks(range(len(imm2))); ax.set_xticklabels(imm2.index, fontsize=9)
ax.set_title("Immune marker scores by subtype")
fig.colorbar(im, ax=ax, shrink=0.8)

ax = fig.add_subplot(gs[1, 1:])
hg = hub.gene.tolist()
Z = np.vstack([(X[idx[g]] - X[idx[g]].mean()) / X[idx[g]].std() for g in hg])
means = pd.DataFrame(Z.T, columns=hg)
means["subtype"] = lab
hm = means.groupby("subtype").mean().T
im = ax.imshow(hm.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
ax.set_yticks(range(len(hg))); ax.set_yticklabels(hg, fontsize=8)
ax.set_xticks(range(len(hm.columns))); ax.set_xticklabels(hm.columns)
ax.set_title("Hub gene expression (mean z) by subtype")
fig.colorbar(im, ax=ax, shrink=0.8)
fig.tight_layout(); fig.savefig(FIG / "fig5_subtypes.png", dpi=200); plt.close(fig)

print("figures saved")