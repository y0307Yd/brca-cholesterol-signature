# -*- coding: utf-8 -*-
"""Build SMR .ma files (SNP A1 A2 freq b se p n) from the BCAC OncoArray
public release (Michailidou et al. 2017) for overall, ER+ and ER- disease,
using the combined OncoArray + iCOGS + GWAS meta columns.

The release's phase3_1kg_id column holds chr:pos:ref:alt variant IDs, so rs
IDs are recovered by matching (chr, pos, allele set) against the 1000 Genomes
EUR PLINK .bim used as the SMR LD reference."""
import gzip
import time
from collections import defaultdict
from pathlib import Path

SRC = r"data\smr\bcac\oncoarray_bcac_public_release_oct17.txt.gz"
BIM = r"data\smr\g1000_eur\g1000_eur.bim"

# Sample sizes (Michailidou 2017 Nature combined meta; EUR)
# overall: 122,977 cases + 105,974 controls; ER+: 69,501 + 105,974;
# ER-: 21,468 + 105,974.
SPECS = {
    "BCAC_Onco_overall.ma": ("bcac_onco_icogs_gwas_", 228951),
    "BCAC_Onco_ERpos.ma": ("bcac_onco_icogs_gwas_erpos_", 175475),
    "BCAC_Onco_ERneg.ma": ("bcac_onco_icogs_gwas_erneg_", 127442),
}


def main():
    bim_map = defaultdict(list)
    with open(BIM, encoding="utf-8") as f:
        for line in f:
            p = line.split()
            bim_map[(int(p[0]), int(p[3]))].append((p[1], p[4], p[5]))
    print("bim position keys:", len(bim_map))

    outs = {name: open(r"data\smr\bcac\\" + name, "w", encoding="utf-8")
            for name in SPECS}
    for name, fh in outs.items():
        fh.write("SNP A1 A2 freq b se p n\n")

    t0 = time.time()
    counts = {name: 0 for name in SPECS}
    skipped = 0
    with gzip.open(SRC, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            v = line.rstrip("\n").split("\t")
            chr_ = v[idx["chr"]].strip()
            pos = v[idx["position_b37"]].strip()
            a0 = v[idx["a0"]].strip()
            a1 = v[idx["a1"]].strip()
            if not (chr_ and pos and a0 and a1):
                skipped += 1
                continue
            key = (int(chr_), int(pos))
            cands = bim_map.get(key, [])
            pal = {a0, a1} in ({"A", "T"}, {"C", "G"})
            matches = []
            for rs, ba1, ba2 in cands:
                if {ba1, ba2} != {a0, a1}:
                    continue
                if pal:
                    # ambiguous strand: require exact allele order
                    if (a0 == ba1 and a1 == ba2) or (a0 == ba2 and a1 == ba1):
                        matches.append((rs, ba1, ba2))
                else:
                    matches.append((rs, ba1, ba2))
            if len(matches) != 1:
                skipped += 1
                continue
            rs = matches[0][0]
            for name, (prefix, n) in SPECS.items():
                eaf, beta, se, p = (v[idx[prefix + "eaf_controls"]],
                                    v[idx[prefix + "beta"]],
                                    v[idx[prefix + "se"]],
                                    v[idx[prefix + "P1df"]])
                if any(x in ("", "NA", ".", "null") for x in (eaf, beta, se, p)):
                    continue
                outs[name].write(
                    f"{rs} {a1} {a0} {eaf} {beta} {se} {p} {n}\n")
                counts[name] += 1
    for name, fh in outs.items():
        fh.close()
    print("rows:", counts, "skipped:", skipped, "sec:", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
