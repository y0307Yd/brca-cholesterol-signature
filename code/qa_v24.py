# -*- coding: utf-8 -*-
"""QA checks for v24 (normal-vs-tumour + immunotherapy analyses)."""
import re

import docx

PATH = (r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature"
        r"\Manuscript_Cholesterol_Metabolism_BRCA_v24.docx")

doc = docx.Document(PATH)
paras = [p.text.strip() for p in doc.paragraphs if p.text.strip()]

i = paras.index("Abstract")
abs_text = " ".join(paras[i + 1:i + 5])
print("Abstract words:", len(abs_text.split()))

ri = paras.index("References")
body = " ".join(paras[5:ri])
body = re.sub(r"Figure Legends.*", "", body, flags=re.S)
print("Body words (excl. refs/figures):", len(body.split()))
print("Body paragraphs:", len(paras[5:ri]))

refs = []
for t in paras[ri + 1:]:
    m = re.match(r"^\[(\d+)\]", t)
    if m:
        refs.append(int(m.group(1)))
    elif t.startswith("Figure Legends"):
        break
print("Reference count:", len(refs), "| sequential 1..n:",
      refs == list(range(1, len(refs) + 1)))

cits = re.findall(r"\[(\d+(?:[,-]\d+)*)\]", body)
order = []
for c in cits:
    for part in c.split(","):
        part = part.strip()
        if "-" in part:
            a, b = map(int, part.split("-"))
            order.extend(range(a, b + 1))
        elif part:
            order.append(int(part))
seen = []
for x in order:
    if x not in seen:
        seen.append(x)
print("First-appearance order OK:", seen == list(range(1, 54)), "|",
      "max ref used:", max(seen), "| missing:",
      sorted(set(range(1, 54)) - set(seen)))

caps = [int(re.match(r"Figure (\d+)\.", t).group(1))
        for t in paras if re.match(r"Figure \d+\.", t)]
print("Figure captions 1..16 complete:", caps == list(range(1, 17)))

supps = sorted({int(x) for x in re.findall(r"Supplementary (?:Table|Figure) S(\d+)", body)},
               key=int)
print("Supplementary items mentioned:", supps)

# restored v15 detail checks
restored = {
    "2.4 WGCNA detail": "median absolute deviation" in body,
    "2.5 z-score + CV grid": "z-scored with TCGA statistics" in body
    and "25-point grid" in body,
    "2.8 scRNA resolutions": "29 minor and 58 subset" in body,
    "2.8 SMR params": "MAF > 0.01" in body and "0.00042" in body,
    "2.8 conditional formula": "(z - r*z_cond)/sqrt(1 - r^2)" in body,
    "3.2 slope/M3105": "slope = -1.64" in body and "M3105" in body,
    "3.2 module details": "15 of 55 cholesterol genes" in body,
    "3.3 Cox details": "HR 0.89 per SD" in body and "4.2e-05" in body,
    "3.4 PAC details": "k = 3 PAC 0.113" in body,
    "3.4 METABRIC ER comp": "C1 15.6%" in body,
    "3.5 DEG counts": "28,931" in body and "15,975" in body,
    "3.5 enrichment details": "hsa04110" in body,
    "3.6 scRNA numbers": "mean 0.45, 46.9%" in body,
    "3.7 SMR details": "P = 4.6e-4" in body and "rs6677385" in body,
    "3.8 GTEx breast details": "24,290 probes" in body and "PON1" in body,
    "3.9 coloc details": "F = 93.1" in body and "155,556,971" in body,
    "3.10 FDPS stage": "rho = 0.04, P = 0.19" in body,
    "3.11 GTEx whole-blood detail": "rs12091730: expected z = -1.27" in body,
    "3.12 Cramer V": "Cramer's V = 0.560" in body,
    "3.13 GSE21653 survival": "HR per SD = 1.27" in body,
    "3.14 SPP1-CD44": "SPP1-CD44" in body,
    "3.15 methylation detail": "mean beta 0.018" in body,
    "4.1 five layers": "Each layer of the framework was completed" in body,
    "4.2 hub list": "DHCR7, DHCR24, NSDHL" in body,
    "4.5 limitations full": "WGCNA-guided prioritization" in body,
}
for k, v in restored.items():
    print(("OK  " if v else "MISS"), k)

# v19/v20 additions still present
added = {
    "GSE20711": "GSE20711, n = 88" in body,
    "CIBERSORT/ESTIMATE": "CIBERSORT" in body and "ESTIMATE" in body,
    "bootstrap": "73.7%" in body,
    "checkpoint": "CD274" in body,
    "HPA": "Human Protein Atlas" in body,
    "pancancer": "17 queried cancer types" in body,
    "coloc bonus": "PP.H4 = 0.035" in body,
    "4.7 comparison": "4.7 Comparison with published" in body,
    "TRIPOD": "TRIPOD-AI checklist" in body,
    "METABRIC checkpoints": "all seven immune-checkpoint genes" in body,
    "GSE20711 subtypes": "15.8%, 58.1%, 79.2% and 15.4%" in body,
    "Fig S4": "Supplementary Figure S4" in body,
    "Fig S5": "Supplementary Figure S5" in body,
    "Fig S6": "Supplementary Figure S6" in body,
    "Table S17": "Supplementary Table S17" in body,
    "Table S18": "Supplementary Table S18" in body,
    "no legacy dup": "80%%" not in "\n".join(paras),
    "pathway panel results": "fatty-acid oxidation" in body,
    "Fig S7": "Supplementary Figure S7" in body,
    "4.3 corrected": "does not segregate" in body,
    "old 4.3 removed": "cholesterol-synthesis-high subtypes are luminal-like" not in body,
    "3.16 normal vs tumour": "3.16 Cholesterol genes are remodelled" in body,
    "3.17 immunotherapy": "3.17 Exploratory immunotherapy-response" in body,
    "GSE15852 result": "GSE15852, 43 pairs" in body,
    "GSE91061 result": "GSE91061" in body and "CD8A P = 0.096" in body,
    "Fig S9/S10": "Supplementary Figure S9" in body
    and "Supplementary Figure S10" in body,
    "Table S19/S20": "Supplementary Table S19" in body
    and "Supplementary Table S20" in body,
    "S1-S20/S1-S10": "S1-S20" in "\n".join(paras)
    and "S1-S10" in "\n".join(paras),
}
for k, v in added.items():
    print(("OK  " if v else "MISS"), "[new]", k)

allok = (refs == list(range(1, len(refs) + 1))
         and seen == list(range(1, 54))
         and caps == list(range(1, 17))
         and all(restored.values()) and all(added.values()))
print("\nQA", "PASSED" if allok else "FAILED")
if not allok:
    raise SystemExit(1)
