#!/usr/bin/env python3
"""Integrate native gene candidates and optional historical Trinotate report.
FAIR: one row per original superTranscript, preserve raw one-to-many protein
annotations in a separate TSV, never infer identities by stripping suffixes.
Native exon overlap and Swiss-Prot similarity are complementary, not proof
of orthology. Tan's historical Trinotate v3.1.1 batch explicitly loads
D. melanogaster Swiss-Prot BLASTX/BLASTP, Pfam and TMHMM; SignalP/RNAmmer
loads are commented. Actual report contents and provenance must be audited.
No functional hits are invented or silently replaced by native gene names.
Manuscript Methods describe both full Swiss-Prot and D. melanogaster subset;
keep reports separately labelled and independently audit exact ID matches.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import sys

def read(path):
    with open(path,encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f,delimiter="\t")

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def write(path,fields,rows):
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison",type=Path,required=True)
    p.add_argument("--prepared",type=Path,required=True)
    p.add_argument("--trinotate",type=Path,help="Legacy single report")
    p.add_argument("--trinotate-full",type=Path)
    p.add_argument("--trinotate-dmel",type=Path)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    lut=a.prepared/"superTranscript_ID_lookup.tsv"
    source=list(read(lut))
    ids={r["original_id"] for r in source}
    if len(ids)!=len(source):raise ValueError("Duplicate source IDs")
    paths=[lut,Path(__file__)]
    evidence={}
    for sp in ("D_simulans","D_melanogaster"):
        path=a.comparison/sp/"gene_summary.tsv"
        paths.append(path)
        rr=list(read(path))
        evidence[sp]={r["original_id"]:r for r in rr}
        if len(rr)!=len(ids) or set(evidence[sp])!=ids:
            raise ValueError("Native summary identifier mismatch: "+sp)
    a.out.mkdir(parents=True,exist_ok=False)
    if a.trinotate and (a.trinotate_full or a.trinotate_dmel):
        p.error("Use either legacy --trinotate or the named database reports")
    reports={"trinotate":a.trinotate} if a.trinotate else {"trinotate_full":a.trinotate_full,"trinotate_dmel":a.trinotate_dmel}
    report_data={}
    for label,report in reports.items():
        hits=collections.defaultdict(list)
        unmatched=[]
        report_fields=[]
        if report:
            paths.append(report)
            with report.open(encoding="utf-8-sig",newline="") as f:
                reader=csv.DictReader(f,delimiter="\t")
                report_fields=reader.fieldnames or []
                canonical={x.lstrip("#").strip().lower():x for x in report_fields}
                gene=canonical.get("gene_id")
                tx=canonical.get("transcript_id")
                if not gene or not tx:
                    raise ValueError("Expected Trinotate gene_id and transcript_id headers; found "+repr(report_fields))
                for r in reader:
                    if None in r or any(v is None for v in r.values()):
                        raise ValueError("Malformed Trinotate row "+str(reader.line_num))
                    matches={r[c] for c in (gene,tx) if r[c] in ids}
                    if len(matches)!=1:
                        unmatched.append(dict(r,join_status="ambiguous" if matches else "unmatched"))
                        continue
                    name=next(iter(matches))
                    hits[name].append(dict(r,original_id=name))
            if "original_id" in report_fields or "join_status" in report_fields:
                raise ValueError("Reserved report header collision")
            write(a.out/f"{label}_matched_rows.tsv",["original_id"]+report_fields,(r for rr in hits.values() for r in rr))
            write(a.out/f"{label}_unmatched_rows.tsv",["join_status"]+report_fields,unmatched)
        report_data[label]=(report,hits,unmatched,report_fields)
    master=[]
    for r in source:
        name=r["original_id"]
        x={k:r[k] for k in ("original_id","projection_id","length","block_count")}
        for sp,rr in evidence.items():
            x.update({sp+"_"+k:v for k,v in rr[name].items() if k!="original_id"})
        for label,(report,hits,unmatched,report_fields) in report_data.items():
            x[label+"_status"]="matched" if hits[name] else "no_matching_report_row" if report else "not_supplied"
            x[label+"_rows"]=len(hits[name]) if report else "NA"
            for field in report_fields:
                vals=list(dict.fromkeys(y[field] for y in hits[name] if y[field] not in ("",".","NA")))
                x[label+"_"+field.lstrip("#")]=json.dumps(vals,ensure_ascii=False) if vals else ""
        master.append(x)
    write(a.out/"master_annotation.tsv",list(master[0]),master)
    summary={"source_sequences":len(ids),"master_rows":len(master),"reports":{}}
    for label,(report,hits,unmatched,fields) in report_data.items():
        summary["reports"][label]={"supplied":bool(report),"matched_sequences":sum(bool(v) for v in hits.values()),"matched_rows":sum(map(len,hits.values())),"unmatched_or_ambiguous_rows":len(unmatched)}
    (a.out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (a.out/"provenance.json").write_text(json.dumps({"command":sys.argv,"python":sys.version,"inputs":{str(x.resolve()):sha(x) for x in paths}},indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)

if __name__=="__main__":
    main()
