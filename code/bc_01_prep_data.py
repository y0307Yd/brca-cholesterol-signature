# -*- coding: utf-8 -*-
"""Prep TCGA/METABRIC inputs for the cholesterol-metabolism signature pipeline."""
import re, json
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(r".\work")
OUT = Path(r".\outputs\chol_metab_signature")
ZHE = Path(r"Path(__file__).resolve().parent.parent")
OUT.mkdir(parents=True, exist_ok=True)

# ---------------- TCGA ----------------
counts = pd.read_csv(WORK / "my_deseq2_input" / "my_counts.csv", index_col=0)
ann = pd.read_csv(WORK / "my_deseq2_input" / "my_annotation.csv", index_col=0)
assert counts.shape[1] == 952 and list(counts.columns) == list(ann.index)
genes = counts.index.tolist()
X = counts.to_numpy(dtype=np.float64)          # genes x samples

lib = X.sum(axis=0)
cpm = X / lib * 1e6
Xlog = np.log2(cpm + 1.0).astype(np.float32)

# clinical traits
clin = pd.read_csv(ZHE / "data" / "processed" / "final_analysis" / "tcga_pfi_clinical_final.csv")
clin = clin[clin["has_valid_pfi"]].copy()
clin["event"] = clin["PFI"].astype(int)
clin["time"] = pd.to_numeric(clin["PFI.time"], errors="coerce").astype(float)
sm = {"I": 1, "II": 2, "III": 3, "IV": 4}
def sn(s):
    if not isinstance(s, str): return np.nan
    m = re.search(r"(IV|III|II|I)", s.upper())
    return sm.get(m.group(1), np.nan) if m else np.nan
clin["stage"] = clin["ajcc_pathologic_tumor_stage"].map(sn)
clin["age"] = pd.to_numeric(clin["age_at_initial_pathologic_diagnosis"], errors="coerce")
clin = clin.set_index("patient_id")

# receptor status from bcr XML
rec = {}
for x in (ZHE / "data" / "raw" / "tcga_brca" / "clinical").glob("*/nationwidechildrens.org_clinical.TCGA-*.xml"):
    m = re.search(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x.name)
    if not m: continue
    t = x.read_text(encoding="utf-8", errors="replace")
    def fld(n):
        mm = re.search(r"<(?:\w+:)?%s[^>]*>([^<]+)</(?:\w+:)?%s>" % (n, n), t)
        return mm.group(1).strip() if mm else None
    rec[m.group(0)] = (fld("breast_carcinoma_estrogen_receptor_status"),
                       fld("breast_carcinoma_progesterone_receptor_status"),
                       fld("lab_proc_her2_neu_immunohistochemistry_receptor_status"))

def enc(v):
    if v == "Positive": return 1.0
    if v == "Negative": return 0.0
    return np.nan
def enc_h(v):
    if v == "Positive": return 1.0
    if v == "Negative": return 0.0
    if v == "Equivocal": return 0.5
    return np.nan

traits = pd.DataFrame({
    "sample_id": list(ann.index),
    "patient_id": ann["patient"].values,
    "event": ann["event"].values.astype(int),
})
t2 = traits["patient_id"].map(clin["time"])
traits["time"] = t2.values
traits["stage"] = traits["patient_id"].map(clin["stage"]).values
traits["age"] = traits["patient_id"].map(clin["age"]).values
traits["ER"] = traits["patient_id"].map(lambda p: enc(rec[p][0]) if p in rec else np.nan).values
traits["PR"] = traits["patient_id"].map(lambda p: enc(rec[p][1]) if p in rec else np.nan).values
traits["HER2"] = traits["patient_id"].map(lambda p: enc_h(rec[p][2]) if p in rec else np.nan).values
assert traits["event"].sum() == 122
assert traits["time"].notna().all()

np.save(OUT / "tcga_Xlog.npy", Xlog)
with open(OUT / "tcga_genes.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(genes))
traits.to_csv(OUT / "tcga_traits.csv", index=False, encoding="utf-8-sig")

# ---------------- METABRIC ----------------
d = ZHE / "data" / "processed" / "full_omics"
mb = np.load(d / "metabric_expr.npz", allow_pickle=True)
mbX = mb["X"].astype(np.float64)               # genes x samples
mb_genes = [l.strip() for l in open(d / "metabric_genes.txt", encoding="utf-8") if l.strip()]
mb_rfs = pd.read_csv(d / "metabric_rfs.csv")
idx = pd.read_csv(ZHE / "data" / "processed" / "metabric" / "metabric_patient_analysis_index.csv")
idx = idx.set_index("SAMPLE_ID")
mb_traits = pd.DataFrame({"sample_id": mb_rfs["SAMPLE_ID"], "event": mb_rfs["event"].astype(int)})
mb_traits["time"] = (mb_rfs["time"] * 30.44).values
mb_traits["ER"] = idx.loc[mb_traits["sample_id"], "ER_STATUS"].map(
    lambda v: 1.0 if isinstance(v, str) and "Positive" in v else (0.0 if isinstance(v, str) and "Negative" in v else np.nan)).values
mb_traits["stage"] = pd.to_numeric(idx.loc[mb_traits["sample_id"], "TUMOR_STAGE"], errors="coerce").values
mb_traits.loc[mb_traits["stage"] == 0, "stage"] = np.nan
assert len(mb_traits) == 1979 and mbX.shape[1] == 1979
np.save(OUT / "mb_X.npy", mbX)
with open(OUT / "mb_genes.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(mb_genes))
mb_traits.to_csv(OUT / "mb_traits.csv", index=False, encoding="utf-8-sig")

# shared genes
shared = sorted(set(genes) & set(mb_genes))
with open(OUT / "shared_genes.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(shared))

summary = {
    "tcga_samples": 952, "tcga_events": int(traits["event"].sum()),
    "tcga_genes": len(genes), "mb_samples": 1979, "mb_events": int(mb_traits["event"].sum()),
    "mb_genes": len(mb_genes), "shared_genes": len(shared),
    "tcga_er_known": int(traits["ER"].notna().sum()), "mb_er_known": int(mb_traits["ER"].notna().sum()),
}
json.dump(summary, open(OUT / "prep_summary.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print(json.dumps(summary, ensure_ascii=False))