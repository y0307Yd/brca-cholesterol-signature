import gzip, time
src=r"data\smr\bcac\GCST010098.h.tsv.gz"; dst=r"data\smr\bcac\BCAC_overall.ma"
t0=time.time(); n=0; skipped=0
with gzip.open(src,"rt") as f, open(dst,"w") as o:
    header=f.readline().rstrip("\n").split("\t")
    idx={c:i for i,c in enumerate(header)}
    o.write("SNP A1 A2 freq b se p n\n")
    for line in f:
        v=line.rstrip("\n").split("\t")
        rs=v[idx["rsid"]].strip()
        if not rs or rs=="NA" or rs==".":
            skipped+=1; continue
        a1=v[idx["effect_allele"]]; a2=v[idx["other_allele"]]
        if not a1 or not a2: skipped+=1; continue
        o.write(f"{rs} {a1} {a2} {v[idx['effect_allele_frequency']]} {v[idx['beta']]} {v[idx['standard_error']]} {v[idx['p_value']]} 247173\n")
        n+=1
print("written", n, "rows, skipped", skipped, "in", round(time.time()-t0,1), "s")