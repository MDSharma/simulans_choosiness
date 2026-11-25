#!/usr/bin/env python3
"""Direct projection of Necklace graph blocks through existing primary BAMs.
FAIR / scientific rationale (M. D. Sharma, 2026-11-25):
Necklace graph blocks partition a composite superTranscript, not native exons.
Original PSL-to-GTF represented alignment blocks but lost source block identity.
Liftoff explicit hierarchy fixed extraction (165 -> 15319); at d=2 it
recovered 96807/96012 distinct graph block IDs in simulans/melanogaster.
Here use existing splice-aware Minimap2 CIGARs, without remapping. Preserve
all original IDs, including 165 differing source transcript labels.
Only M/= /X consume both query and reference; I/S/H do not project, D/N
break genomic pieces. Reverse BAM orientation is transformed to original
FASTA coordinates. Primary alignments alone determine the main tracks;
secondary/supplementary counts are retained as ambiguity evidence, not merged
into a potentially chimeric gene. A mapped base need not match exactly.
Coordinates: internal and TSV *_0 intervals are 0-based half-open, GFF3
1-based closed. Output graph_block fragments are not validated native exons.
FAIR: deterministic source block IDs matching prepared GFF3, reversible
lookup, checksums, command and version provenance, explicit unmapped rows,
original files untouched. References: doi:10.1186/s13059-017-1284-1;
doi:10.1093/gigascience/giy045; doi:10.1093/bioinformatics/bty191;
doi:10.1093/bioinformatics/btaa1016; SAM specification at
https://samtools.github.io/hts-specs/SAMv1.pdf.
"""
import argparse
import bisect
import collections
import csv
import datetime
import hashlib
import json
from pathlib import Path
import sys
from urllib.parse import quote
import pysam

def log(s):
    print(datetime.datetime.now().isoformat(timespec="seconds"), s, flush=True)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()

def attrs(s):
    return dict(x.split("=", 1) for x in s.strip(";").split(";") if "=" in x)

def writer(f, fields):
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(fields)
    return w

def pieces(cigar, start, qlen, reverse):
    """Return aligned query/reference runs; hard clips count in original q."""
    q, r = 0, start
    runs = []
    for op, n in cigar:
        if op in (0, 7, 8):
            runs.append((qlen-q-n, qlen-q, r, r+n) if reverse else (q, q+n, r, r+n))
            q += n
            r += n
        elif op in (1, 4, 5):
            q += n
        elif op in (2, 3):
            r += n
        elif op != 6:
            raise ValueError("Unsupported CIGAR operation: " + str(op))
    if q != qlen:
        raise ValueError(f"CIGAR query length {q} != source {qlen}")
    return sorted(runs)

def project(intervals, runs, reverse):
    # Sweep source intervals against aligned query runs. Coalesce only when
    # both query AND reference are adjacent, so insertions/deletions stay visible.
    ends = [x[1] for x in runs]
    for start, end, bid in intervals:
        parts = []
        for i in range(bisect.bisect_right(ends, start), len(runs)):
            a,b,c,d = runs[i]
            if a >= end:
                break
            x,y = max(a,start), min(b,end)
            if x >= y:
                continue
            u,v = (d-(y-a), d-(x-a)) if reverse else (c+x-a,c+y-a)
            if parts:
                p,z,j,k = parts[-1]
                if z == x and ((k == u and not reverse) or (v == j and reverse)):
                    parts[-1] = (p,y,min(j,u),max(k,v))
                    continue
            parts.append((x,y,u,v))
        yield start,end,bid,parts

def self_test():
    # 2S 4M 1I 3M 5N 2M 2D 2M 1H, length 15.
    c = [(4,2),(0,4),(1,1),(0,3),(3,5),(0,2),(2,2),(0,2),(5,1)]
    assert pieces(c,100,15,False) == [(2,6,100,104),(7,10,104,107),(10,12,112,114),(12,14,116,118)]
    assert pieces(c,100,15,True) == [(1,3,116,118),(3,5,112,114),(5,8,104,107),(9,13,100,104)]
    assert list(project([(4,11,"b")],pieces(c,100,15,False),False))[0][3] == [(4,6,102,104),(7,10,104,107),(10,11,112,113)]
    assert list(project([(4,11,"b")],pieces(c,100,15,True),True))[0][3] == [(4,5,112,113),(5,8,104,107),(9,11,102,104)]
    assert list(project([(0,4,"b")],pieces([(7,2),(8,2)],0,4,False),False))[0][3] == [(0,4,0,4)]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mapping-run",type=Path,default=Path("SuperTranscript_mapping_20261125_215819"))
    p.add_argument("--prepared",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--species",nargs="+",choices=["D_simulans","D_melanogaster"],default=["D_simulans","D_melanogaster"])
    args=p.parse_args()
    self_test()
    lut=args.prepared/"superTranscript_ID_lookup.tsv"
    gff=args.prepared/"SuperDuper.blocks.gff3"
    with lut.open() as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    source={r["original_id"]:r for r in rows}
    by_sid={r["projection_id"]:r["original_id"] for r in rows}
    if len(source)!=len(rows) or len(by_sid)!=len(rows):
        raise ValueError("Duplicate source IDs")
    blocks=collections.defaultdict(list)
    with gff.open() as f:
        for l in f:
            if l.startswith("#") or not l.strip(): continue
            a=l.rstrip("\n").split("\t")
            if a[2]=="exon":
                blocks[by_sid[a[0]]].append((int(a[3])-1,int(a[4]),attrs(a[8])["ID"]))
    seen=set()
    for name,row in source.items():
        pos=0
        blocks[name].sort()
        for a,b,bid in blocks[name]:
            if a!=pos or b<=a or bid in seen: raise ValueError("Invalid source block "+bid)
            seen.add(bid)
            pos=b
        if pos!=int(row["length"]) or len(blocks[name])!=int(row["block_count"]):
            raise ValueError("Source partition mismatch "+name)
    # Verify exact sequence identity between prepared FASTA and BAM query source.
    with pysam.FastaFile(str(args.mapping_run/"tmp/query.fa")) as orig, pysam.FastaFile(str(args.prepared/"SuperDuper.standardised.fasta")) as prep:
        if set(orig.references)!=set(source) or set(prep.references)!=set(by_sid):
            raise ValueError("Source FASTA identifier mismatch")
        for name,row in source.items():
            if orig.fetch(name)!=prep.fetch(row["projection_id"]):
                raise ValueError("Source sequence mismatch "+name)
    args.out.mkdir(parents=True,exist_ok=False)
    provenance={"command":sys.argv,"python":sys.version,"pysam":pysam.__version__,
                "inputs":{str(x.resolve()):sha(x) for x in [lut,gff,args.mapping_run/"tmp/query.fa",args.prepared/"SuperDuper.standardised.fasta",Path(__file__)]}}
    with (args.out/"summary.tsv").open("w") as sf:
        sw=writer(sf,["species","queries","primary_mapped","source_blocks","full_blocks","partial_blocks","unmapped_blocks","source_bases","projected_query_bases","fragments"])
        for species in args.species:
            log(f"{species}: scanning existing BAM")
            bam=args.mapping_run/"alignments"/f"{species}_all.bam"
            provenance["inputs"][str(bam.resolve())]=sha(bam)
            dest=args.out/species
            dest.mkdir()
            primary={}
            counts=collections.defaultdict(lambda:[0,0])
            with pysam.AlignmentFile(str(bam),"rb") as f:
                for a in f.fetch(until_eof=True):
                    name=a.query_name
                    if name not in source: raise ValueError("Unknown BAM query "+name)
                    if a.is_unmapped: continue
                    if a.is_secondary: counts[name][0]+=1
                    elif a.is_supplementary: counts[name][1]+=1
                    else:
                        if name in primary: raise ValueError("Multiple primary alignments "+name)
                        primary[name]=(a.reference_name,a.reference_start,a.cigartuples,a.is_reverse,a.mapping_quality)
            total=collections.Counter()
            features=[]
            with (dest/"blocks.tsv").open("w") as bf,(dest/"fragments.tsv").open("w") as ff,(dest/"queries.tsv").open("w") as qf:
                bw=writer(bf,["original_id","block_id","query_start_0","query_end_0","mapped_query_bases","status","fragments","contig","strand","mapq"])
                fw=writer(ff,["original_id","block_id","fragment_id","query_start_0","query_end_0","contig","ref_start_0","ref_end_0","strand"])
                qw=writer(qf,["original_id","length","primary_mapped","mapped_query_bases","full_blocks","partial_blocks","unmapped_blocks","secondary_alignments","supplementary_alignments","contig","strand","mapq"])
                for name,row in source.items():
                    qlen=int(row["length"])
                    if name in primary:
                        contig,start,cigar,rev,mapq=primary[name]
                        runs=pieces(cigar,start,qlen,rev)
                        strand="-" if rev else "+"
                    else:
                        contig,strand,mapq="NA","NA","NA"
                        runs,rev=[],False
                    qc=collections.Counter()
                    for a,b,bid,parts in project(blocks[name],runs,rev):
                        mapped=sum(y-x for x,y,u,v in parts)
                        status="full" if mapped==b-a else "partial" if mapped else "unmapped"
                        qc[status]+=1
                        qc["bases"]+=mapped
                        total["fragments"]+=len(parts)
                        bw.writerow([name,bid,a,b,mapped,status,len(parts),contig,strand,mapq])
                        for i,(x,y,u,v) in enumerate(parts,1):
                            fid=f"{bid}.part{i}"
                            fw.writerow([name,bid,fid,x,y,contig,u,v,strand])
                            at=f"ID={fid};Name={quote(name,safe=':._-')};original_id={quote(name,safe=':._-')};block_id={bid};query_start={x+1};query_end={y};block_status={status}"
                            features.append((contig,u,v,f"{contig}\tNecklace_CIGAR\tregion\t{u+1}\t{v}\t.\t{strand}\t.\t{at}\n"))
                    total.update(qc)
                    qw.writerow([name,qlen,int(name in primary),qc["bases"],qc["full"],qc["partial"],qc["unmapped"],*counts[name],contig,strand,mapq])
            with (dest/"Necklace_CIGAR_blocks.gff3").open("w") as f:
                f.write("##gff-version 3\n")
                for _,_,_,line in sorted(features): f.write(line)
            sw.writerow([species,len(source),len(primary),len(seen),total["full"],total["partial"],total["unmapped"],sum(int(r["length"]) for r in rows),total["bases"],total["fragments"]])
            sf.flush()
            log(f"{species}: {len(primary)} primary mapped; blocks {total}")
    (args.out/"provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    log("Complete: "+str(args.out))

if __name__=="__main__":
    main()
