#!/usr/bin/env bash
# SuperTranscript genome mapping; M. D. Sharma, revised pipeline.
set -Eeuo pipefail
WORKDIR=$(realpath "${WORKDIR:-$PWD}")
cd "$WORKDIR"
SUPERTRANSCRIPTS=$(realpath "${SUPERTRANSCRIPTS:-superTranscriptome/SuperDuper.fasta}")
SUPERTRANSCRIPT_GTF=$(realpath "${SUPERTRANSCRIPT_GTF:-superTranscriptome/SuperDuper.gtf}")
SIM_GENOME=$(realpath "${SIM_GENOME:-Simulans_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fa}")
MEL_GENOME=$(realpath "${MEL_GENOME:-Melanogaster_genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fa}")
SIM_GTF=$(realpath "${SIM_GTF:-Simulans_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.gtf}")
MEL_GTF=$(realpath "${MEL_GTF:-Melanogaster_genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf}")
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
LIFT_DISTANCES=${LIFT_DISTANCES:-'2 50'}
LIFT_MM2_OPTIONS=${LIFT_MM2_OPTIONS:-}
RUN_LIFTOFF=${RUN_LIFTOFF:-1}
DE_LIST=${DE_LIST:-}
OUTROOT=$(realpath -m "${OUTROOT:-SuperTranscript_mapping_$(date +%Y%m%d_%H%M%S)}")
[[ $THREADS =~ ^[0-9]+$ ]] && (( THREADS >= 4 )) || { echo 'Need at least 4 CPUs' >&2; exit 1; }
[[ ! -e $OUTROOT ]] || { echo "Output already exists: $OUTROOT (choose a new OUTROOT)" >&2; exit 1; }
mkdir -p "$OUTROOT"/{logs,alignments,bed,bigwig,liftoff,qc,annotation,jbrowse_tracks,tmp}
exec > >(tee -a "$OUTROOT/logs/pipeline.log") 2>&1
trap 'echo "FAILED at line $LINENO" >&2' ERR
for t in samtools minimap2 bedtools python3 bgzip tabix; do command -v "$t" >/dev/null || { echo "Missing dependency: $t"; exit 1; }; done
if [[ $RUN_LIFTOFF == 1 ]]; then
    command -v liftoff >/dev/null || { echo "Missing dependency: liftoff. Activate the environment containing Liftoff or install it with mamba install -c conda-forge -c bioconda liftoff." >&2; exit 1; }
fi
for f in "$SUPERTRANSCRIPTS" "$SIM_GENOME" "$MEL_GENOME" "$SIM_GTF" "$MEL_GTF"; do [[ -s $f ]] || { echo "Missing $f"; exit 1; }; done
if [[ $RUN_LIFTOFF == 1 ]]; then [[ -s $SUPERTRANSCRIPT_GTF ]]; fi
if [[ -n $DE_LIST ]]; then DE_LIST=$(realpath "$DE_LIST"); [[ -s $DE_LIST ]]; fi
{
 echo "Started $(date -Is)"; echo "THREADS=$THREADS LIFT_DISTANCES=$LIFT_DISTANCES LIFT_MM2_OPTIONS=$LIFT_MM2_OPTIONS"
 minimap2 --version; samtools --version; bedtools --version
 if [[ $RUN_LIFTOFF == 1 ]]; then liftoff -V; fi
 sha256sum "$SUPERTRANSCRIPTS" "$SIM_GENOME" "$MEL_GENOME" "$SIM_GTF" "$MEL_GTF"
 if [[ $RUN_LIFTOFF == 1 ]]; then sha256sum "$SUPERTRANSCRIPT_GTF"; fi
} > "$OUTROOT/logs/provenance.txt"
# Copy inputs used by index-writing tools, to avoid modifying originals.
for pair in "query:$SUPERTRANSCRIPTS" "D_simulans:$SIM_GENOME" "D_melanogaster:$MEL_GENOME"; do
 name=${pair%%:*}; f=${pair#*:}; cp "$f" "$OUTROOT/tmp/$name.fa"; samtools faidx "$OUTROOT/tmp/$name.fa"
done
QUERY="$OUTROOT/tmp/query.fa"
# Validate sequence IDs and bounds; distinguish source vs genomic annotation.
cat > "$OUTROOT/tmp/validate.py" <<'PY'
import sys, collections
fai, annotation = sys.argv[1:]
lengths = {x.split('\t')[0]:int(x.split('\t')[1]) for x in open(fai)}
counts=collections.Counter(); bad=[]
for n,line in enumerate(open(annotation),1):
    if line.startswith('##FASTA'): break
    if line.startswith('#') or not line.strip(): continue
    f=line.rstrip('\n').split('\t')
    try:
        assert len(f)==9 and f[0] in lengths and 1<=int(f[3])<=int(f[4])<=lengths[f[0]]
    except (AssertionError,ValueError):
        bad.append((n,line[:150])); continue
    counts[f[2]]+=1
print(annotation, dict(counts), 'invalid rows:',len(bad))
if bad: print(*bad[:10],sep='\n')
if bad or not counts: sys.exit(1)
PY
if [[ $RUN_LIFTOFF == 1 ]]; then
 python3 "$OUTROOT/tmp/validate.py" "$QUERY.fai" "$SUPERTRANSCRIPT_GTF" || { echo 'Liftoff source annotation must use SuperDuper FASTA coordinates. Set SUPERTRANSCRIPT_GTF to the original Necklace annotation or RUN_LIFTOFF=0.'; exit 1; }
 cp "$SUPERTRANSCRIPT_GTF" "$OUTROOT/tmp/source.${SUPERTRANSCRIPT_GTF##*.}"
 SOURCE="$OUTROOT/tmp/source.${SUPERTRANSCRIPT_GTF##*.}"
fi
if [[ -n $DE_LIST ]]; then
 python3 - "$QUERY.fai" "$DE_LIST" "$OUTROOT/tmp/de.ids" <<'PY'
import sys
fai,src,out=sys.argv[1:]; known={l.split('\t')[0] for l in open(fai)}
ids=list(dict.fromkeys(l.strip() for l in open(src) if l.strip()))
missing=set(ids)-known
if not ids or missing: sys.exit(f'Empty or unknown DE IDs: {sorted(missing)[:20]}')
open(out,'w').write('\n'.join(ids)+'\n')
PY
fi
SORT_THREADS=1
MAP_THREADS=$((THREADS-3))
for species in D_simulans D_melanogaster; do
 echo "Mapping $species"
 genome="$OUTROOT/tmp/$species.fa"
 if [[ $species == D_simulans ]]; then refgtf=$SIM_GTF; else refgtf=$MEL_GTF; fi
 python3 "$OUTROOT/tmp/validate.py" "$genome.fai" "$refgtf"
 bam="$OUTROOT/alignments/${species}_all.bam"
 minimap2 -ax splice:hq -u b --secondary=yes -N 5 --eqx -t "$MAP_THREADS" "$genome" "$QUERY" | samtools sort -@ "$SORT_THREADS" -m 1G -o "$bam"
 samtools quickcheck "$bam"; samtools index "$bam"
 samtools flagstat "$bam" > "$OUTROOT/qc/${species}_flagstat.txt"
 # One row per alignment; explicit query coverage and identity (excluding introns).
 samtools view "$bam" | python3 -c '
import sys,re
lengths={r.split("\t")[0]:int(r.split("\t")[1]) for r in open(sys.argv[1])}
print("query\tflag\tcontig\tstart\tmapq\tquery_length\taligned_query_bases\tquery_coverage\tidentity")
for line in sys.stdin:
    f=line.rstrip().split("\t"); flag=int(f[1]); qlen=lengths[f[0]]
    ops=[(int(n),op) for n,op in re.findall(r"(\d+)([MIDNSHP=X])",f[5])]
    a=sum(n for n,o in ops if o in "MI=X")
    denom=sum(n for n,o in ops if o in "MID=X")
    nm=next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")),None)
    identity=f"{1-nm/denom:.6f}" if denom and nm is not None else "NA"
    print(f[0],flag,f[2],f[3],f[4],qlen,a,f"{a/qlen:.6f}",identity,sep="\t")
' "$QUERY.fai" > "$OUTROOT/qc/${species}_alignments.tsv"

 primary="$OUTROOT/alignments/${species}_primary.bam"
 samtools view -b -F 2308 -o "$primary" "$bam"; samtools index "$primary"
 track="$OUTROOT/jbrowse_tracks/$species"; mkdir -p "$track"
 cp "$genome" "$genome.fai" "$primary" "$primary.bai" "$track/"
 awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && NF==9' "$refgtf" | LC_ALL=C sort -k1,1 -k4,4n | bgzip -c > "$track/reference.gtf.gz"
 tabix -p gff "$track/reference.gtf.gz"
 for mode in all primary; do
  b="$OUTROOT/alignments/${species}_${mode}.bam"
  bedtools bamtobed -bed12 -i "$b" > "$OUTROOT/bed/${species}_${mode}.bed12"
  bedtools bamtobed -split -i "$b" > "$OUTROOT/bed/${species}_${mode}_blocks.bed"
 done
 cp "$OUTROOT/bed/${species}_primary.bed12" "$track/"
 # Raw alignment-block depth, not expression; split at N and D.
 bedtools genomecov -ibam "$primary" -bg -split > "$OUTROOT/bigwig/${species}.bedgraph"
 if command -v bedGraphToBigWig >/dev/null; then
  cut -f1,2 "$genome.fai" > "$OUTROOT/tmp/$species.sizes"
  LC_ALL=C sort -k1,1 -k2,2n "$OUTROOT/bigwig/${species}.bedgraph" > "$OUTROOT/tmp/$species.sorted.bedgraph"
  if [[ -s $OUTROOT/tmp/$species.sorted.bedgraph ]]; then
   bedGraphToBigWig "$OUTROOT/tmp/$species.sorted.bedgraph" "$OUTROOT/tmp/$species.sizes" "$track/${species}_alignment_depth.bw"
  fi
 fi
 # Candidate reference gene overlaps; keep raw rows for isoform-aware aggregation.
 awk -F '\t' '$3=="exon"' "$refgtf" > "$OUTROOT/tmp/$species.exons.gtf"
 bedtools intersect -a "$OUTROOT/bed/${species}_primary_blocks.bed" -b "$OUTROOT/tmp/$species.exons.gtf" -wo > "$OUTROOT/annotation/${species}_exon_overlaps.tsv"
 if [[ -n $DE_LIST ]]; then
  samtools view -b -N "$OUTROOT/tmp/de.ids" -o "$track/${species}_DE.bam" "$primary"
  samtools index "$track/${species}_DE.bam"
  bedtools bamtobed -bed12 -i "$track/${species}_DE.bam" > "$track/${species}_DE.bed12"
 fi
 if [[ $RUN_LIFTOFF == 1 ]]; then
  for d in $LIFT_DISTANCES; do
   [[ $d =~ ^[0-9]+([.][0-9]+)?$ ]] || exit 1
   dest="$OUTROOT/liftoff/${species}_d${d}"; mkdir -p "$dest"
   opts=(); if [[ -n $LIFT_MM2_OPTIONS ]]; then opts+=("-mm2_options=$LIFT_MM2_OPTIONS"); fi
   # Infer parents only for GTF when absent. GFF must have valid hierarchy.
   if [[ $SOURCE == *.gtf ]]; then
    if ! awk -F '\t' '$3=="gene" {found=1} END{exit !found}' "$SOURCE"; then opts+=(-infer_genes); fi
    if ! awk -F '\t' '$3=="transcript" || $3=="mRNA" {found=1} END{exit !found}' "$SOURCE"; then opts+=(-infer_transcripts); fi
   fi
   if (cd "$dest"; liftoff "$genome" "$QUERY" -g "$SOURCE" -o projected.annotation -u unmapped.txt -dir intermediate -p "$THREADS" -d "$d" -a 0.5 -s 0.5 "${opts[@]}") > "$dest/liftoff.log" 2>&1; then
    # Liftoff preserves source format; detect from output rather than suffix.
    if grep -m1 -q '^##gff-version 3' "$dest/projected.annotation"; then ext=gff3; else ext=${SOURCE##*.}; fi
    cp "$dest/projected.annotation" "$dest/projected.$ext"
    if python3 "$OUTROOT/tmp/validate.py" "$genome.fai" "$dest/projected.$ext"; then
     cp "$dest/projected.$ext" "$track/Necklace_d${d}.$ext"
     if [[ $ext == gtf ]] && command -v gffread >/dev/null; then
      gffread -E -F --keep-genes "$dest/projected.$ext" -o "$track/Necklace_d${d}.converted.gff3" 2> "$dest/gffread.log" || { rm -f "$track/Necklace_d${d}.converted.gff3"; echo "GFF conversion failed; native projection retained"; }
     fi
     printf '%s\t%s\t%s\n' "$species" "$d" OK >> "$OUTROOT/qc/liftoff_status.tsv"
    else printf '%s\t%s\t%s\n' "$species" "$d" INVALID >> "$OUTROOT/qc/liftoff_status.tsv"; fi
   else printf '%s\t%s\t%s\n' "$species" "$d" FAILED >> "$OUTROOT/qc/liftoff_status.tsv"; fi
  done
 fi
done
printf 'Finished %s\nOutputs: %s\n' "$(date -Is)" "$OUTROOT"
if [[ $RUN_LIFTOFF == 1 ]] && grep -Eq 'FAILED|INVALID' "$OUTROOT/qc/liftoff_status.tsv"; then echo 'One or more Liftoff trials failed; inspect logs. Mapping outputs retained.'; exit 2; fi
