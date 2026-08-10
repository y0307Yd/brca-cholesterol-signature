# -*- coding: utf-8 -*-
"""Prepare SuSiE/coloc inputs for the FDPS locus.

Reads the per-SNP eQTLGen FDPS cis table, the BCAC overall GWAS (.ma) and
the 1000 Genomes EUR PLINK bfiles; harmonises alleles, computes z-scores
and an LD correlation matrix, and writes binary inputs for R.
"""
import csv
import struct
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(r"outputs\chol_metab_signature\susie_fdps")
EQTL = OUT.parent / "fdps_eqltgen_cis_snps.csv"
GWAS = Path(r"data\smr\bcac\BCAC_overall.ma")
BFILE = Path(r"data\smr\g1000_eur\g1000_eur")


def load_plink_rows(prefix, snp_indices):
    """Load genotypes for selected SNP rows (0-based) from a PLINK1 bed."""
    fam = [l.split() for l in open(str(prefix) + ".fam", encoding="utf-8")]
    bim = [l.split() for l in open(str(prefix) + ".bim", encoding="utf-8")]
    n_ind = len(fam)
    nbytes = (n_ind + 3) // 4
    geno = np.empty((len(snp_indices), n_ind), dtype=np.int8)
    with open(str(prefix) + ".bed", "rb") as f:
        assert f.read(3) == b"\x6c\x1b\x01"
        for out_i, snp_i in enumerate(snp_indices):
            f.seek(3 + snp_i * nbytes)
            b = np.frombuffer(f.read(nbytes), dtype=np.uint8)
            bits = np.unpackbits(b, bitorder="little")[: n_ind * 2]
            b0 = bits[0::2].astype(np.int8)
            b1 = bits[1::2].astype(np.int8)
            g = b0 + b1
            g[(b0 == 0) & (b1 == 1)] = -9
            geno[out_i] = g
    return bim, geno


def main():
    e = pd.read_csv(EQTL, dtype={"SNP": str, "chr": str, "bp": int})
    gwas = {}
    with open(GWAS, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter=" ")
        for r in rd:
            gwas[r["SNP"]] = r
    print("eQTL SNPs:", len(e), "| GWAS rows:", len(gwas))

    bim = [l.split() for l in open(str(BFILE) + ".bim", encoding="utf-8")]
    bim_idx = {b[1]: i for i, b in enumerate(bim)}
    print("1000G SNPs:", len(bim))

    rows = []
    snp_list = []
    bim_hits = []
    ld_flip = []
    for _, r in e.iterrows():
        rs = r["SNP"]
        if rs not in gwas or rs not in bim_idx:
            continue
        bi = bim_idx[rs]
        b_g = gwas[rs]
        # harmonise: effect allele = eQTL A1
        a1e = r["A1"]
        if b_g["A1"] == a1e:
            zg = float(b_g["b"]) / float(b_g["se"])
        elif b_g["A2"] == a1e:
            zg = -float(b_g["b"]) / float(b_g["se"])
        else:
            continue
        ze = float(r["b_eQTL"]) / float(r["se_eQTL"])
        # LD dosage of the effect allele (a1e) in 1000G
        a1g, a2g = bim[bi][4], bim[bi][5]
        if a2g == a1e:
            ld_flip.append(False)
        elif a1g == a1e:
            ld_flip.append(True)
        else:
            continue
        rows.append([rs, r["chr"], int(r["bp"]), a1e, r["A2"],
                     float(r["freq"]), ze, zg])
        snp_list.append(rs)
        bim_hits.append(bim_idx[rs])
    print("overlap SNPs (eQTL ∩ GWAS ∩ 1000G):", len(snp_list))

    # LD matrix from dosages (mean-imputed)
    _, G = load_plink_rows(BFILE, bim_hits)
    G = G.astype(float)
    for i, fl in enumerate(ld_flip):
        if fl:
            G[i] = np.where(G[i] == -9, -9, 2 - G[i])
    G[G == -9] = np.nan
    means = np.nanmean(G, axis=1)
    for i in range(G.shape[0]):
        m = np.isnan(G[i])
        G[i, m] = means[i]
    # standardise and flip to effect-allele dosage
    # G was coded as 0/1/2 of effect allele for matched rows already
    Gs = (G - G.mean(axis=1, keepdims=True)) / G.std(axis=1, ddof=1, keepdims=True)
    R = (Gs @ Gs.T) / (Gs.shape[1] - 1)
    R = np.nan_to_num(R, nan=0.0)
    np.fill_diagonal(R, 1.0)
    print("LD matrix:", R.shape)
    Z = np.array([[r[6], r[7]] for r in rows], dtype=float)

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["SNP", "chr", "bp", "A1", "A2", "freq",
                                "z_eQTL", "z_GWAS"]).to_csv(
        OUT / "snp_table.csv", index=False, encoding="utf-8-sig")
    np.save(OUT / "ld_matrix.npy", R.astype(np.float32))
    np.save(OUT / "z_matrix.npy", Z.astype(np.float32))
    # binary float64 dump for R
    with open(OUT / "ld_matrix.bin", "wb") as f:
        f.write(struct.pack("<ii", *R.shape))
        f.write(R.astype(np.float64).tobytes())
    with open(OUT / "z_matrix.bin", "wb") as f:
        f.write(struct.pack("<ii", *Z.shape))
        f.write(Z.astype(np.float64).tobytes())
    print("saved to", OUT)


if __name__ == "__main__":
    main()
