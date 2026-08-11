# -*- coding: utf-8 -*-
"""Build Supplementary Tables S13-S15 (v19 bonus items)."""
import json
from pathlib import Path

import pandas as pd

OUT = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature")
WORK = Path(r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\work")

# S13: HPA protein evidence + pan-cancer FDPS expression
hpa = json.load(open(WORK / "hpa_fdps.json", encoding="utf-8"))
bc_tcga = hpa["Cancer prognostics - Breast Invasive Carcinoma (TCGA)"]
bc_val = hpa["Cancer prognostics - Breast Invasive Carcinoma (validation)"]
pc = pd.read_csv(WORK / "bonus_pancancer_fdps.csv")
rows = [
    {"item": "UniProt ID", "value": ";".join(hpa.get("Uniprot", []))},
    {"item": "HPA antibody", "value": ";".join(hpa.get("Antibody", []))},
    {"item": "Antibody reliability (IH)", "value": hpa.get("Reliability (IH)", "")},
    {"item": "Subcellular main location", "value": ";".join(hpa.get("Subcellular main location", []))},
    {"item": "Subcellular additional location", "value": ";".join(hpa.get("Subcellular additional location", []))},
    {"item": "Protein tissue specificity", "value": hpa.get("Protein tissue specificity", "")},
    {"item": "RNA tissue specificity", "value": hpa.get("RNA tissue specificity", "")},
    {"item": "BC TCGA protein prognosis", "value": f"unprognostic, P={bc_tcga.get('p_val')}"},
    {"item": "BC validation protein prognosis", "value": f"unprognostic, P={bc_val.get('p_val')}"},
    {"item": "Biological process", "value": "; ".join(hpa.get("Biological process", [])[:8])},
]
hpa_rows = pd.DataFrame(rows)
pc_rows = pd.DataFrame({
    "item": "pan-cancer expression: " + pc["cancer"],
    "value": "detected (n = " + pc["n"].astype(str) + ")",
})
pc_rows.loc[len(pc_rows)] = ["pan-cancer summary",
                             f"detected in all {len(pc)} cancer types with available RNA-seq"]
combined = pd.concat([hpa_rows, pc_rows], ignore_index=True)
combined.to_csv(OUT / "Supplementary_Table_S13_hpa_pancancer_FDPS.csv",
                index=False, encoding="utf-8-sig")

# S14: CIBERSORT/ESTIMATE/immune-checkpoint by subtype
cib = pd.read_csv(WORK / "bonus_cibersort_by_subtype.csv")
est = pd.read_csv(WORK / "bonus_estimate_by_subtype.csv")
chk = pd.read_csv(WORK / "bonus_checkpoint_by_subtype.csv")
cib = cib.rename(columns={"cell": "cell_type"})
est = est.rename(columns={"score": "score_type"})
chk = chk.rename(columns={"gene": "checkpoint_gene",
                          "highest_subtype": "highest_subtype",
                          "lowest_subtype": "lowest_subtype"})
with pd.ExcelWriter(OUT / "Supplementary_Table_S14_immune_deconvolution.xlsx") as xw:
    cib.to_excel(xw, sheet_name="CIBERSORT_fractions", index=False)
    est.to_excel(xw, sheet_name="ESTIMATE_scores", index=False)
    chk.to_excel(xw, sheet_name="Immune_checkpoints", index=False)

# S15: coloc bonus results (LIPA/FAXDC2/SREBF1)
coloc = pd.read_csv(OUT / "coloc_bonus_summary.csv")
coloc.to_csv(OUT / "Supplementary_Table_S15_coloc_bonus.csv",
             index=False, encoding="utf-8-sig")

print("S13:", combined.shape)
print("S14 sheets written; CIBERSORT rows:", len(cib),
      "ESTIMATE rows:", len(est), "checkpoint rows:", len(chk))
print("S15:", coloc.shape)
