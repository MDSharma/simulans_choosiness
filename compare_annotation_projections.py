#!/usr/bin/env python3
"""Compare direct Necklace block projections with Liftoff and native genes.
FAIR rationale: native exon overlap, Liftoff projection and translated
protein similarity answer different questions. Neither locus overlap nor a
best protein hit proves orthology. Use union of genomic aligned bases,
not repeated exon rows across isoforms. Report both strands; alignment
orientation may differ from biological strand for composite sequences.
Use original Necklace IDs for downstream DE and Trinotate. No remapping.
Coordinates internal 0-based half-open; inputs/outputs checksummed.
Historical benchmark: Liftoff d2 projected 14503/14315 clusters and
96807/96012 distinct graph-block IDs (simulans/melanogaster); all 15319
source clusters extracted after explicit hierarchy fix. Compare block IDs,
not merely counts; Liftoff spans may contain introns or deleted bases.
References: Davidson et al. doi:10.1186/s13059-017-1284-1;
Shumate & Salzberg doi:10.1093/bioinformatics/btaa1016;
Li doi:10.1093/bioinformatics/bty191.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

def rows(path):
    with open(path) as f:
        yield from csv.DictReader(f,delimiter="\t")

def table(path,header,data):
    with open(path,"w") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(header)
        w.writerows(data)

def merge(xs):
    out=[]
    for a,b in sorted(xs):
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(b,out[-1][1]))
        else: out.append((a,b))
    return out

def size(xs):
    return sum(b-a for a,b in merge(xs))

def overlap(a,b):
    a,b=merge(a),merge(b)
    i=j=total=0
    while i<len(a) and j<len(b):
        x,y=a[i]; u,v=b[j]
        total+=max(0,min(y,v)-max(x,u))
        if y<=v: i+=1
        else: j+=1
    return total

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--direct",type=Path,required=True)
    p.add_argument("--mapping-run",type=Path,default=Path("SuperTranscript_mapping_20261125_215819"))
    p.add_argument("--sim-liftoff",type=Path,required=True,help="simulans benchmark run root")
    p.add_argument("--mel-liftoff",type=Path,default=Path("Necklace_annotation_projection_20261125_225413"))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    provenance={"command":sys.argv,"python":sys.version,"bedtools":subprocess.check_output(["bedtools","--version"],text=True).strip(),"inputs":{}}
    def record(path):
        provenance["inputs"][str(path.resolve())]=sha(path)
    record(Path(__file__))
    summaries=[]
    recovered={}
    for sp,run in [("D_simulans",a.sim_liftoff),("D_melanogaster",a.mel_liftoff)]:
        print(sp+": comparing block coordinates and native gene candidates",flush=True)
        d=a.direct/sp
        out=a.out/sp
        out.mkdir()
        blockfile=d/"blocks.tsv"; fragfile=d/"fragments.tsv"; queryfile=d/"queries.tsv"
        liftfile=run/"liftoff"/sp/"blocks_d2/projected.gff3"
        liftgenes=run/"liftoff"/sp/"blocks_d2/per_gene_projection.tsv"
        lutfile=run/"prepared_inputs/superTranscript_ID_lookup.tsv"
        ref=a.mapping_run/"jbrowse_tracks"/sp/"reference.gtf.gz"
        for path in [blockfile,fragfile,queryfile,liftfile,liftgenes,lutfile,ref]: record(path)
        blocks={r["block_id"]:r for r in rows(blockfile)}
        lut={r["original_id"]:r for r in rows(lutfile)}
        queries={r["original_id"]:r for r in rows(queryfile)}
        # Verify each block ID corresponds to the same source cluster and index.
        for bid,r in blocks.items():
            sid=lut[r["original_id"]]["projection_id"]
            if not bid.startswith("graph_block_"+sid+"_"):
                raise ValueError("Incompatible prepared identifier lookup: "+bid)
        direct=collections.defaultdict(lambda:collections.defaultdict(list))
        lift=collections.defaultdict(lambda:collections.defaultdict(list))
        with (out/"aligned_fragments.bed").open("w") as f:
            for r in rows(fragfile):
                key=(r["contig"],r["strand"])
                x,y=int(r["ref_start_0"]),int(r["ref_end_0"])
                direct[r["block_id"]][key].append((x,y))
                f.write("\t".join([key[0],str(x),str(y),r["original_id"],"0",key[1]])+"\n")
        with liftfile.open() as f:
            for line in f:
                if line.startswith("#") or not line.strip():continue
                r=line.rstrip("\n").split("\t")
                if r[2]!="exon":continue
                attr=dict(x.split("=",1) for x in r[8].strip(";").split(";") if "=" in x)
                bid=unquote(attr["ID"])
                if bid not in blocks: raise ValueError("Unknown Liftoff block "+bid)
                lift[bid][(r[0],r[6])].append((int(r[3])-1,int(r[4])))
        recovered[sp]=set(lift)
        counts=collections.Counter()
        data=[]
        for bid,r in blocks.items():
            ds,ls=direct.get(bid,{}),lift.get(bid,{})
            db=sum(size(v) for v in ds.values()); lb=sum(size(v) for v in ls.values())
            ov=sum(overlap(v,ls[k]) for k,v in ds.items() if k in ls)
            same_contig=bool({k[0] for k in ds}&{k[0] for k in ls})
            state="both" if ds and ls else "direct_only" if ds else "liftoff_only" if ls else "neither"
            counts[state]+=1
            counts["same_strand_overlap"]+=ov>0
            data.append([r["original_id"],bid,r["status"],state,int(same_contig),int(bool(set(ds)&set(ls))),db,lb,ov,ov/db if db else "NA",ov/lb if lb else "NA"])
        table(out/"block_comparison.tsv",["original_id","block_id","direct_status","recovery","same_contig","same_contig_strand","direct_ref_bases","liftoff_ref_span_bases","overlap_same_strand_bases","fraction_direct_overlap","fraction_liftoff_overlap"],data)
        summaries.append([sp,len(blocks),counts["both"],counts["direct_only"],counts["liftoff_only"],counts["neither"],counts["same_strand_overlap"]])
        # Extract native exon intervals and metadata, with a stable gene ID.
        meta={}
        def parse(s):
            return dict(re.findall(r'(\S+)\s+"([^"]*)"',s))
        with gzip.open(ref,"rt") as f,(out/"native_exons.bed").open("w") as g:
            for line in f:
                if line.startswith("#") or not line.strip():continue
                r=line.rstrip("\n").split("\t")
                if r[2] not in ("gene","exon"):continue
                at=parse(r[8]); gid=at.get("gene_id")
                if not gid: raise ValueError("Missing native gene_id")
                key=(gid,r[0],r[6])
                m=meta.setdefault(key,{})
                for field in ("gene","gene_name","gene_biotype","gene_type","db_xref","description"):
                    if at.get(field):m.setdefault(field,at[field])
                if r[2]=="exon":
                    g.write("\t".join([r[0],str(int(r[3])-1),r[4],gid,"0",r[6]])+"\n")
        # Bedtools intersect emits overlap for each isoform exon; union per
        # query/gene/locus removes duplicate counting of shared exon bases.
        candidates=collections.defaultdict(list)
        cmd=["bedtools","intersect","-a",str(out/"aligned_fragments.bed"),"-b",str(out/"native_exons.bed"),"-wo"]
        with subprocess.Popen(cmd,stdout=subprocess.PIPE,text=True) as proc:
            for l in proc.stdout:
                r=l.rstrip("\n").split("\t")
                key=(r[3],r[9],r[0],r[5],r[11])
                candidates[key].append((max(int(r[1]),int(r[7])),min(int(r[2]),int(r[8]))))
            if proc.wait():raise RuntimeError("bedtools intersect failed")
        by_query=collections.defaultdict(list)
        for (name,gid,contig,qs,gs),intervals in candidates.items():
            n=size(intervals)
            denom=int(queries[name]["mapped_query_bases"])
            m=meta.get((gid,contig,gs),{})
            by_query[name].append([name,gid,contig,qs,gs,int(qs==gs),n,n/denom if denom else 0,m.get("gene",m.get("gene_name","")),m.get("gene_biotype",m.get("gene_type","")),m.get("db_xref",""),m.get("description","")])
        header=["original_id","reference_gene_id","contig","alignment_strand","gene_strand","same_strand","union_exon_overlap_bp","fraction_projected_query_bases","gene_symbol","gene_biotype","db_xref","description"]
        ordered=[]
        gene_summary=[]
        lf={r["original_id"]:r for r in rows(liftgenes)}
        for name,q in queries.items():
            c=sorted(by_query[name],key=lambda r:(-r[6],-r[5],r[1]))
            ordered.extend(c)
            same=[r for r in c if r[5]]
            # Explicitly retain ties, no forced one-to-one gene assignment.
            best=[r for r in same if r[6]==same[0][6]] if same else []
            gene_summary.append([name,q["primary_mapped"],q["mapped_query_bases"],q["mapq"],q["secondary_alignments"],q["supplementary_alignments"],len(c),len(same),";".join(r[1] for r in best),";".join(r[8] for r in best),best[0][6] if best else 0,best[0][7] if best else 0,lf[name]["projected"],lf[name]["primary_span_overlap"]])
        table(out/"reference_gene_candidates.tsv",header,ordered)
        table(out/"gene_summary.tsv",["original_id","primary_mapped","mapped_query_bases","mapq","secondary_alignments","supplementary_alignments","reference_candidates_all","reference_candidates_same_strand","top_same_strand_gene_ids","top_same_strand_gene_symbols","top_overlap_bp","top_overlap_fraction","liftoff_projected","liftoff_primary_span_overlap"],gene_summary)
        print(sp+": "+str(counts),flush=True)
    table(a.out/"block_comparison_summary.tsv",["species","source_blocks","both","direct_only","liftoff_only","neither","same_strand_overlap"],summaries)
    s,m=recovered["D_simulans"],recovered["D_melanogaster"]
    table(a.out/"liftoff_shared_blocks.tsv",["category","count"],[["both",len(s&m)],["simulans_only",len(s-m)],["melanogaster_only",len(m-s)],["neither",len(blocks)-len(s|m)]])
    (a.out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    print("Complete: "+str(a.out),flush=True)

if __name__=="__main__":
    main()
