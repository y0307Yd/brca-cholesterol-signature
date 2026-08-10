# -*- coding: utf-8 -*-
"""Step 2: LASSO + SVM-RFE hub genes, ER-status classifier, external validation."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
from lifelines import CoxPHFitter

OUT = Path(r".\outputs\chol_metab_signature")
RNG = 20260804
np.random.seed(RNG)

X = np.load(OUT / "tcga_Xlog.npy")
genes = [l.strip() for l in open(OUT / "tcga_genes.txt", encoding="utf-8") if l.strip()]
tr = pd.read_csv(OUT / "tcga_traits.csv")
mbX = np.load(OUT / "mb_X.npy")
mb_genes = [l.strip() for l in open(OUT / "mb_genes.txt", encoding="utf-8") if l.strip()]
mb = pd.read_csv(OUT / "mb_traits.csv")

shared = set(open(OUT / "shared_genes.txt", encoding="utf-8").read().splitlines())
chol = pd.read_csv(OUT / "cholesterol_genes.csv")["symbol"].tolist()
cand = sorted(set(chol) & shared)
idx = {g: i for i, g in enumerate(genes)}
midx = {g: i for i, g in enumerate(mb_genes)}
print("candidate genes (cholesterol x shared):", len(cand))

Xc = np.vstack([X[idx[g]] for g in cand]).T          # samples x genes
Mc = np.vstack([mbX[midx[g]] for g in cand]).T
# z-score with TCGA statistics
mu = Xc.mean(axis=0); sd = Xc.std(axis=0); sd[sd == 0] = 1.0
Xz = (Xc - mu) / sd
Mz = (Mc - mu) / sd

mask = tr["ER"].notna().values
y = tr["ER"].values[mask].astype(int)
Xtr = Xz[mask]

# ---- LASSO C selection by 5-fold CV ----
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG)
C_grid = np.logspace(-3, 1, 25)
cv_scores = []
for C in C_grid:
    aucs = []
    for tr_i, va_i in cv.split(Xtr, y):
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=2000, random_state=RNG)
        clf.fit(Xtr[tr_i], y[tr_i])
        aucs.append(roc_auc_score(y[va_i], clf.predict_proba(Xtr[va_i])[:, 1]))
    cv_scores.append(np.mean(aucs))
C_best = float(C_grid[int(np.argmax(cv_scores))])
print("best C:", C_best, "cv AUC:", round(max(cv_scores), 4))

clf = LogisticRegression(penalty="l1", solver="liblinear", C=C_best, max_iter=2000, random_state=RNG)
clf.fit(Xtr, y)
lasso_genes = [cand[i] for i, c in enumerate(clf.coef_[0]) if abs(c) > 1e-6]
print("LASSO selected:", len(lasso_genes), lasso_genes)

# ---- SVM-RFE ----
nsel = min(len(lasso_genes), 12)
if nsel < 3:
    nsel = min(3, len(cand))
rfe = RFE(estimator=SVC(kernel="linear", C=1.0, class_weight="balanced", random_state=RNG),
          n_features_to_select=nsel, step=1)
rfe.fit(Xtr, y)
rfe_genes = [cand[i] for i in range(len(cand)) if rfe.support_[i]]
print("SVM-RFE selected:", len(rfe_genes), rfe_genes)

hub = sorted(set(lasso_genes) & set(rfe_genes))
if len(hub) < 2:
    hub = sorted(set(lasso_genes) | set(rfe_genes))[:8]
print("hub genes:", hub)

# ---- final model on hub genes ----
hx = {g: cand.index(g) for g in hub}
Xh = Xz[:, [hx[g] for g in hub]]
Mh = Mz[:, [hx[g] for g in hub]]
clf2 = LogisticRegression(penalty="l1", solver="liblinear", C=C_best, max_iter=2000, random_state=RNG)
clf2.fit(Xh[mask], y)
coef = dict(zip(hub, clf2.coef_[0]))

# TCGA internal CV AUC (hub genes, C fixed)
aucs = []
for tr_i, va_i in cv.split(Xh[mask], y):
    m = LogisticRegression(penalty="l1", solver="liblinear", C=C_best, max_iter=2000, random_state=RNG)
    m.fit(Xh[mask][tr_i], y[tr_i])
    aucs.append(roc_auc_score(y[va_i], m.predict_proba(Xh[mask][va_i])[:, 1]))
cv_auc_mean = float(np.mean(aucs)); cv_auc_sd = float(np.std(aucs))
full_auc = roc_auc_score(y, clf2.predict_proba(Xh[mask])[:, 1])

# METABRIC external
my = mb["ER"].astype(int).values
p_mb = clf2.predict_proba(Mh)[:, 1]
ext_auc = roc_auc_score(my, p_mb)
boot = []
rng = np.random.default_rng(RNG)
for _ in range(2000):
    ix = rng.integers(0, len(my), len(my))
    if len(np.unique(my[ix])) == 2:
        boot.append(roc_auc_score(my[ix], p_mb[ix]))
ext_ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
print("TCGA CV AUC %.3f ± %.3f | full %.3f | METABRIC %.3f (95%% CI %.3f-%.3f)" %
      (cv_auc_mean, cv_auc_sd, full_auc, ext_auc, ext_ci[0], ext_ci[1]))

# signature score & prognosis (secondary, honest)
sig_tcga = Xh @ clf2.coef_[0]
sig_mb = Mh @ clf2.coef_[0]
cox_rows = []
for name, t, e, s in [("TCGA", tr["time"].values, tr["event"].values.astype(int), sig_tcga),
                      ("METABRIC", mb["time"].values, mb["event"].values.astype(int), sig_mb)]:
    d = pd.DataFrame({"time": t, "event": e, "sig": s})
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(d, duration_col="time", event_col="event", show_progress=False)
    r = cph.summary.iloc[0]
    cox_rows.append({"cohort": name, "HR": float(r["exp(coef)"]), "p": float(r["p"]),
                     "ci_low": float(r["exp(coef) lower 95%"]), "ci_high": float(r["exp(coef) upper 95%"])})
cox = pd.DataFrame(cox_rows)
print(cox.to_string(index=False))

# save
pd.DataFrame({"gene": hub, "coef": [coef[g] for g in hub]}).to_csv(
    OUT / "hub_genes.csv", index=False, encoding="utf-8-sig")
np.save(OUT / "signature_tcga.npy", sig_tcga)
np.save(OUT / "signature_mb.npy", sig_mb)
fpr, tpr, _ = roc_curve(my, p_mb)
pd.DataFrame({"fpr": fpr, "tpr": tpr}).to_csv(OUT / "roc_metabric.csv", index=False)
res = {"C_best": C_best, "lasso_genes": lasso_genes, "rfe_genes": rfe_genes,
       "hub_genes": hub, "cv_auc_mean": cv_auc_mean, "cv_auc_sd": cv_auc_sd,
       "full_auc_tcga": full_auc, "external_auc_metabric": ext_auc,
       "external_auc_ci95": ext_ci, "cox": cox_rows}
json.dump(res, open(OUT / "ml_summary.json", "w", encoding="utf-8"), indent=1)
print("saved step2 outputs")