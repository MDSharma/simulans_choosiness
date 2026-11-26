#!/usr/bin/env python3
"""Project recovered blocks.gtf independently of historical count ordering.
Reuse audited CIGAR functions; do not infer limma IDs or biological exons.
Source row IDs are stable for this exact GTF, not statistical feature IDs.
"""
import argparse, csv, collections, hashlib, json
from pathlib import Path
from urllib.parse import quote
import pysam
from project_necklace_blocks import pieces, project, self_test
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--mapping-run",type=Path,required=True)
p.add_argument("--blocks",type=Path,required=True)
p.add_argument("--out",type=Path,required=True)
a=p.parse_args()
self_test()
blocks=collections.defaultdict(list)
with pysam.FastaFile(str(a.mapping_run/"tmp/query.fa")) as f:
    lengths=dict(zip(f.references,f.lengths))
with a.blocks.open() as f:
    for n,line in enumerate(f,1):
        if line.startswith("#") or not line.strip(): continue
        r=line.rstrip("\n").split("\t")
        if len(r)!=9: raise ValueError(f"Invalid GTF row {n}")
        if r[2]!="exon": continue
        name=r[0]; x,y=int(r[3])-1,int(r[4])
        if not 0<=x<y<=lengths[name]: raise ValueError(f"Invalid coordinates row {n}")
        blocks[name].append((x,y,f"counted_row_{n:07d}"))
if not blocks: raise ValueError("No exon rows")
for v in blocks.values(): v.sort()
a.out.mkdir(parents=True,exist_ok=False)
with (a.out/"source_blocks.tsv").open("w") as f:
    w=csv.writer(f,delimiter="\t",lineterminator="\n")
    w.writerow(["original_id","block_id","start_1based","end_1based"])
    for name,v in blocks.items():
        for x,y,b in v: w.writerow([name,b,x+1,y])
summary={}
for sp in ("D_simulans","D_melanogaster"):
    primary={}
    with pysam.AlignmentFile(str(a.mapping_run/"alignments"/f"{sp}_all.bam"),"rb") as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            if r.query_name in primary: raise ValueError("Duplicate primary")
            primary[r.query_name]=(r.reference_name,r.reference_start,r.cigartuples,r.is_reverse)
    counts=collections.Counter()
    with (a.out/f"{sp}.gff3").open("w") as g,(a.out/f"{sp}.tsv").open("w") as t:
        g.write("##gff-version 3\n")
        w=csv.writer(t,delimiter="\t",lineterminator="\n")
        w.writerow(["original_id","block_id","source_start_1based","source_end_1based","projected_bases","status","contig","orientation","fragments"])
        for name,v in blocks.items():
            contig,start,cigar,rev=primary.get(name,("NA",0,[],False))
            runs=pieces(cigar,start,lengths[name],rev) if cigar else []
            strand="-" if rev else "+"
            for x,y,b,parts in project(v,runs,rev):
                mapped=sum(j-i for i,j,u,z in parts)
                status="full" if mapped==y-x else "partial" if mapped else "unmapped"
                counts[status]+=1
                w.writerow([name,b,x+1,y,mapped,status,contig,strand if cigar else ".",len(parts)])
                for k,(i,j,u,z) in enumerate(parts,1):
                    at=f"ID={b}.part{k};Name={quote(name,safe=':._-')};original_id={quote(name,safe=':._-')};block_id={b};source_start={x+1};source_end={y};query_start={i+1};query_end={j};block_status={status}"
                    g.write(f"{contig}\tNecklace_counted\tregion\t{u+1}\t{z}\t.\t{strand}\t.\t{at}\n")
    summary[sp]=dict(counts)
    print(sp,dict(counts),flush=True)
(a.out/"summary.json").write_text(json.dumps({"source_gtf_sha256":hashlib.sha256(a.blocks.read_bytes()).hexdigest(),"status":summary,"note":"Structural projection only; statistical row correspondence unverified"},indent=2)+"\n")
