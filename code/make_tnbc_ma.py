import gzip, time
src=r"data\smr\bcac\GCST010100.h.tsv.gz"; dst=r"data\smr\bcac\BCAC_TNBC.ma"
t0=time.time(); n=0; skip=0
with gzip.open(src,"rt") as f, open(dst,"w") as o:
    h=f.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(h)}
    o.write("SNP A1 A2 freq b se p n\n")
    for line in f:
        v=line.rstrip("\n").split("\t")
        rs=v[ix["rsid"]].strip()
        a1=v[ix["effect_allele"]]; a2=v[ix["other_allele"]]
        if not rs or rs in ("NA",".") or not a1 or not a2:
            skip+=1; continue
        o.write(f"{rs} {a1} {a2} {v[ix['effect_allele_frequency']]} {v[ix['beta']]} {v[ix['standard_error']]} {v[ix['p_value']]} 118987\n")
        n+=1
print("rows", n, "skipped", skip, "sec", round(time.time()-t0,1))