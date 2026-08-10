"""Merge HEIDI P-values from the dedicated HEIDI run into Table S2."""
from pathlib import Path

import pandas as pd

OUT = Path(r"outputs\chol_metab_signature")


def load(path, sep):
    return pd.read_csv(path, sep=sep, encoding="utf-8-sig")


pathway = load(OUT / "smr_chol_pathway_results.csv", ",")
heidi = load(Path(r"data\smr\results\smr_heidi.smr"), "\t")
heidi_map = dict(zip(heidi["probeID"], heidi["p_HEIDI"]))

orig = dict(zip(pathway["probeID"], pathway["p_HEIDI"]))


def merge(p):
    v = heidi_map.get(p)
    if v and v not in ("NA", ""):
        return v
    return orig.get(p, "NA")


pathway["p_HEIDI"] = pathway["probeID"].map(merge)
pathway.to_csv(OUT / "Supplementary_Table_S2_SMR_pathway_results.csv",
               index=False, encoding="utf-8-sig")
fdps = pathway[pathway["Gene"] == "ENSG00000160752"]
print(fdps[["probeID", "p_SMR", "p_HEIDI"]].to_string(index=False))
