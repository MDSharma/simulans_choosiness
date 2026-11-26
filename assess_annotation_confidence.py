#!/usr/bin/env python3
"""Assess genomic support for Necklace gene associations, not orthology.
FAIR rationale: PSL conversion loses original block identity; explicit Liftoff
hierarchies resolved 165/15319 extraction bias; direct CIGAR projection preserves
block identities but primary placements and query orientation are not proof of
transcriptional strand. Retain all candidates and exact original IDs. Thresholds
are exploratory, not calibrated probabilities. Native exon overlap is union
across isoforms; denominator is projected query bases, not total source length.
Report LOC tokens as evidence only, not assumed gene IDs. No remapping.
References: Davidson et al. 2017 doi:10.1186/s13059-017-1284-1;
Davidson & Oshlack 2018 doi:10.1093/gigascience/giy045;
Li 2018 doi:10.1093/bioinformatics/bty191;
Shumate & Salzberg 2021 doi:10.1093/bioinformatics/btaa1016.
"""
import argparse
import collections
import csv
import hashlib
import json
import pathlib
import re
import sys

def read(p):
    with p.open() as f:
        yield from csv.DictReader(f, delimiter="\t")

def write(p, fields, rows):
    with p.open("w") as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison",type=pathlib.Path,default=pathlib.Path("Necklace_projection_comparison"))
    p.add_argument("--mapping-run",type=pathlib.Path,default=pathlib.Path("SuperTranscript_mapping_20261125_215819"))
    p.add_argument("--master",type=pathlib.Path,default=pathlib.Path("Necklace_integrated_annotation_native/master_annotation.tsv"))
    p.add_argument("--out",type=pathlib.Path,required=True)
    a=p.parse_args()
    master=list(read(a.master))
    ids={r["original_id"] for r in master}
    if len(ids)!=len(master): raise ValueError("Duplicate master IDs")
    a.out.mkdir(parents=True,exist_ok=False)
    inputs=[a.master,pathlib.Path(__file__)]
    summary=[]
    for sp in ("D_simulans","D_melanogaster"):
        print("Assessing "+sp,flush=True)
        qpath=a.mapping_run/"qc"/(sp+"_alignments.tsv")
        cpath=a.comparison/sp/"reference_gene_candidates.tsv"
        inputs.extend([qpath,cpath])
        primary={}
        for r in read(qpath):
            flag=int(r["flag"])
            if not flag & (4|256|2048):
                if r["query"] in primary: raise ValueError("Duplicate primary")
                primary[r["query"]]=r
        candidates=collections.defaultdict(list)
        for r in read(cpath): candidates[r["original_id"]].append(r)
        counts=collections.Counter()
        for r in master:
            name=r["original_id"]
            pre=sp+"_"
            q=primary.get(name)
            c=candidates[name]
            same=[x for x in c if x["same_strand"]=="1"]
            frac=float(r[pre+"top_overlap_fraction"])
            tied=len(r[pre+"top_same_strand_gene_ids"].split(";"))>1
            cov=float(q["query_coverage"]) if q else 0
            mq=int(q["mapq"]) if q else 0
            alt=int(r[pre+"secondary_alignments"])
            supp=int(r[pre+"supplementary_alignments"])
            if not q: category="G_unmapped"
            elif not c: category="F_no_exon_overlap"
            elif not same: category="E_opposite_orientation_only"
            elif tied or len(same)>1 or supp or alt: category="D_ambiguous"
            elif mq>=20 and cov>=0.8 and frac>=0.8: category="A_strong"
            elif mq>=20 and cov>=0.5 and frac>=0.5: category="B_moderate"
            else: category="C_limited"
            counts[category]+=1
            tokens=set(re.findall(r"LOC\d+",name))
            observed=set()
            for x in c:
                observed.update(re.findall(r"LOC\d+"," ".join([x["reference_gene_id"],x["gene_symbol"],x["db_xref"]])))
            r.update({pre+"support_category":category,pre+"primary_contig":q["contig"] if q else "NA",
                      pre+"primary_start_1based":q["start"] if q else "NA",
                      pre+"primary_query_coverage":q["query_coverage"] if q else "NA",
                      pre+"primary_identity":q["identity"] if q else "NA",
                      pre+"source_LOC_tokens":";".join(sorted(tokens)),
                      pre+"overlapping_LOC_tokens":";".join(sorted(tokens&observed)),
                      pre+"multiple_candidate_flag":int(len(same)>1),pre+"tied_top_flag":int(tied)})
        summary.extend({"species":sp,"category":k,"count":v} for k,v in sorted(counts.items()))
    write(a.out/"master_with_confidence.tsv",list(master[0]),master)
    write(a.out/"confidence_summary.tsv",["species","category","count"],summary)
    (a.out/"provenance.json").write_text(json.dumps({"command":sys.argv,"python":sys.version,"inputs":{str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs}},indent=2)+"\n")
    print("Complete: "+str(a.out),flush=True)

if __name__=="__main__": main()
