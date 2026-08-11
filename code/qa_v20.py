# -*- coding: utf-8 -*-
"""QA checks for v20 manuscript."""
import re

import docx

PATH = (r"C:\Users\Y\Documents\Codex\2026-08-06\new-chat\outputs\chol_metab_signature"
        r"\Manuscript_Cholesterol_Metabolism_BRCA_v20.docx")

doc = docx.Document(PATH)
paras = [p.text.strip() for p in doc.paragraphs if p.text.strip()]

# abstract word count
i = paras.index("Abstract")
abs_text = " ".join(paras[i + 1:i + 5])
print("Abstract words:", len(abs_text.split()))

# body: everything between Abstract keywords and References
ri = paras.index("References")
body = " ".join(paras[5:ri])
body = re.sub(r"Figure Legends.*", "", body, flags=re.S)
print("Body words (excl. refs/figures):", len(body.split()))

# citation order
refs = []
for t in paras[ri + 1:]:
    m = re.match(r"^\[(\d+)\]", t)
    if m:
        refs.append(int(m.group(1)))
    elif t.startswith("Figure Legends"):
        break
print("Reference count:", len(refs), "| sequential 1..n:",
      refs == list(range(1, len(refs) + 1)))

# first appearance order of citations in body text
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

# figure mentions
figs = sorted({int(x) for x in re.findall(r"Figure (\d+)", body)})
print("Figures mentioned:", figs, "| 1..16 complete:", figs == list(range(1, 17)))

# supplementary mentions
supps = sorted({x for x in re.findall(r"Supplementary (?:Table|Figure) S(\d+)", body)},
               key=int)
print("Supplementary items mentioned:", supps)

# new content presence
checks = {
    "LASSO bootstrap": "bootstrap resamples" in body,
    "CIBERSORT": "CIBERSORT" in body,
    "ESTIMATE": "ESTIMATE" in body,
    "checkpoint": "Immune-checkpoint genes" in body,
    "GEO subtype validation": "GSE21653 (n = 252) and " in body,
    "GSE20711 results": "GSE20711, n = 88" in body,
    "Fackler ref": "[43]" in body and "Fackler" in body,
    "HPA": "Human Protein Atlas" in body,
    "pan-cancer": "17 queried cancer types" in body,
    "coloc bonus": "PP.H4 = 0.035" in body,
    "comparison 4.7": "4.7 Comparison with published" in body,
    "TRIPOD": "TRIPOD-AI checklist" in body,
    "HPA ref [42]": "[42]" in body,
    "new refs": all(f"[{n}]" in "\n".join(paras) for n in [48, 49, 50, 51, 52, 53]),
    "old limitation removed": "not CIBERSORT/ESTIMATE [40,41]" not in body,
}
for k, v in checks.items():
    print(("OK  " if v else "MISS"), k)

if not all(checks.values()):
    raise SystemExit("QA failed")
print("\nQA PASSED")
