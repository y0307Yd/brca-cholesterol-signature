import csv
from pathlib import Path

import pandas as pd

RES = Path(r"data\smr\results")
OUT = Path(r"outputs\chol_metab_signature")


def load(tag):
    rows = []
    with open(RES / f"smr_chol_onco_{tag}.smr", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["analysis"] = tag
            rows.append(r)
    return rows


frames = []
for tag in ["ERpos", "ERneg"]:
    frames.append(pd.DataFrame(load(tag)))
df = pd.concat(frames, ignore_index=True)
df = df.sort_values(["analysis", "p_SMR"])
df.to_csv(OUT / "Supplementary_Table_S10_ER_stratified_SMR.csv",
          index=False, encoding="utf-8-sig")
print("rows:", len(df), "| FDPS rows:")
print(df[df["Gene"] == "ENSG00000160752"].to_string(index=False))
