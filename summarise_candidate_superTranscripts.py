#!/usr/bin/env python3
"""Extract manuscript candidates by exact original ID; preserve all evidence.
FAIR: explicit candidate list, no guessing from gene symbols or fuzzy joins.
DESeq2 and limma diffSplice results are distinct; do not infer exon skipping
from differential block usage. Candidate TSV must contain original_id and
may include any statistics or labels, including repeated rows for blocks.
Source graph blocks (142593) differ from counted blocks (121757 reported).
Retain raw result rows and native candidate tables; no remapping.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

def read(p):
    with p.open() as f:
        r=csv.DictReader(f,delimiter="\t")
        if not r.fieldnames or "original_id" not in r.fieldnames:
            raise ValueError("Expected original_id column in "+str(p))
        return r.fieldnames,list(r)

def write(p,h,rr):
    with p.open("w") as f:
        w=csv.DictWriter(f,fieldnames=h,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rr)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master",type=Path,default=Path("Necklace_annotation_confidence/master_with_confidence.tsv"))
    p.add_argument("--candidates",type=Path,required=True)
    p.add_argument("--comparison",type=Path,default=Path("Necklace_projection_comparison"))
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    mh,m=read(a.master)
    ch,c=read(a.candidates)
    ids={r["original_id"] for r in c}
    known={r["original_id"] for r in m}
    if not ids or ids-known: raise ValueError("Empty/unknown candidate IDs: "+str(sorted(ids-known)[:20]))
    a.out.mkdir(parents=True,exist_ok=False)
    write(a.out/"candidate_master.tsv",mh,(r for r in m if r["original_id"] in ids))
    write(a.out/"candidate_input.tsv",ch,c)
    inputs=[a.master,a.candidates,Path(__file__)]
    for sp in ("D_simulans","D_melanogaster"):
        path=a.comparison/sp/"reference_gene_candidates.tsv"
        h,rr=read(path)
        write(a.out/(sp+"_reference_gene_candidates.tsv"),h,(r for r in rr if r["original_id"] in ids))
        inputs.append(path)
    (a.out/"candidate_ids.txt").write_text("\n".join(sorted(ids))+"\n")
    (a.out/"provenance.json").write_text(json.dumps({"command":sys.argv,"python":sys.version,"inputs":{str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs}},indent=2)+"\n")
    print(f"Retained {len(ids)} distinct candidates from {len(c)} input rows",flush=True)
if __name__=="__main__": main()
