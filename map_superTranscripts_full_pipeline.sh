#!/usr/bin/env bash
# SuperTranscript genome mapping; M. D. Sharma, revised pipeline.
# ----------------------------------------------------------------------
# SCIENTIFIC RATIONALE AND FAIR-ORIENTED PIPELINE DESIGN
# Author: M. D. Sharma | ORCID: 0000-0002-9957-3153
# Repository: https://github.com/MDSharma/simulans_choosiness
#
# Objective: locate Necklace SuperDuper sequences against D. simulans and
# D. melanogaster, preserve alignment evidence and compare native annotations.
# The original BLAT PSL -> GTF converter correctly converted alignment
# coordinates but did not project original Necklace features, distinguish
# multiple loci, or establish genuine exon boundaries. Direct splice-aware
# Minimap2 gives BAM/BED12 and QC; Liftoff is an annotation-aware complement,
# not a proven replacement or a guarantee of orthology or native isoforms.
#
# Necklace/Lace: SuperDuper.fasta is a linear union of sequence per inferred
# gene cluster. SuperDuper.gtf represents whole sequences for counting;
# SuperDuper.gff contains graph-derived blocks. Neither is a conventional
# genomic gene annotation. Splice-junction-derived blocks are distinct.
# Historical run (2026-11-25, SuperTranscript_mapping_20261125_215819):
# 15319 sequences; 15064 and 14940 primary mappings to simulans and
# melanogaster, respectively. Median primary query coverage 99.92% / 99.02%;
# median identity 99.51% / 93.37%. Liftoff extracted only 165 sequences
# (163 REL-), independently of distance factor or target. This was an
# extraction/hierarchy problem, not evidence of widespread mapping failure.
# Exact cause remains unconfirmed. Do not interpret the historical Liftoff
# status OK (coordinate validation only) as complete annotation recovery.
# Use prepare_necklace_annotation.py and evaluate_liftoff_projection.sh
# for explicit hierarchy, reversible IDs, extraction audits and controlled
# whole-vs-block projection comparisons without repeating direct mapping.
# The operational code below is deliberately retained for reproducibility.
#
# Findable: named species, stage-specific outputs, timestamped progress.
# Accessible: open formats FASTA, BAM/BAI, BED12, GTF/GFF3, TSV, BigWig;
# access permissions, licensing and archival deposition are separate decisions.
# Interoperable: original Necklace IDs for downstream expression/Trinotate;
# reference-coordinate tracks kept separately by species.
# Reusable: tool versions, input SHA256, retained intermediates, exact
# parameters, original inputs preserved. Raw alignment depth is not RNA
# expression; reference exon overlaps are candidate evidence, not definitive
# gene assignments. FAIR-oriented design is not a certification.
#
# References:
# Davidson et al. 2017 doi:10.1186/s13059-017-1284-1
# Davidson & Oshlack 2018 doi:10.1093/gigascience/giy045
# Shumate & Salzberg 2021 doi:10.1093/bioinformatics/btaa1016
# Li 2018 doi:10.1093/bioinformatics/bty191
# https://github.com/Oshlack/necklace/wiki/Output
# https://github.com/Oshlack/Lace/wiki/Usage-Documentation
# ----------------------------------------------------------------------
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
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
trap 'rc=$?; log "FAILED at line $LINENO (exit $rc): $BASH_COMMAND" >&2' ERR
log "Starting superTranscript mapping | CPUs=$THREADS | output=$OUTROOT"
log "Checking software dependencies and input files"
for t in samtools minimap2 bedtools python3 bgzip tabix; do command -v "$t" >/dev/null || { echo "Missing dependency: $t"; exit 1; }; done
if [[ $RUN_LIFTOFF == 1 ]]; then
    command -v liftoff >/dev/null || { echo "Missing dependency: liftoff. Activate the environment containing Liftoff or install it with mamba install -c conda-forge -c bioconda liftoff." >&2; exit 1; }
fi
for f in "$SUPERTRANSCRIPTS" "$SIM_GENOME" "$MEL_GENOME" "$SIM_GTF" "$MEL_GTF"; do [[ -s $f ]] || { echo "Missing $f"; exit 1; }; done
if [[ $RUN_LIFTOFF == 1 ]]; then [[ -s $SUPERTRANSCRIPT_GTF ]]; fi
if [[ -n $DE_LIST ]]; then DE_LIST=$(realpath "$DE_LIST"); [[ -s $DE_LIST ]]; fi
log "Recording tool versions and input SHA256 checksums (may take a moment)"
{
 echo "Started $(date -Is)"; echo "THREADS=$THREADS LIFT_DISTANCES=$LIFT_DISTANCES LIFT_MM2_OPTIONS=$LIFT_MM2_OPTIONS"
 minimap2 --version; samtools --version; bedtools --version
 if [[ $RUN_LIFTOFF == 1 ]]; then liftoff -V; fi
 sha256sum "$SUPERTRANSCRIPTS" "$SIM_GENOME" "$MEL_GENOME" "$SIM_GTF" "$MEL_GTF"
 if [[ $RUN_LIFTOFF == 1 ]]; then sha256sum "$SUPERTRANSCRIPT_GTF"; fi
} > "$OUTROOT/logs/provenance.txt"
log "Provenance recorded in $OUTROOT/logs/provenance.txt"
log "Preparing indexed working copies of FASTA inputs"
# Copy inputs used by index-writing tools, to avoid modifying originals.
for pair in "query:$SUPERTRANSCRIPTS" "D_simulans:$SIM_GENOME" "D_melanogaster:$MEL_GENOME"; do
 name=${pair%%:*}; f=${pair#*:}; log "Copying and indexing $name: $f"; cp "$f" "$OUTROOT/tmp/$name.fa"; samtools faidx "$OUTROOT/tmp/$name.fa"
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
 log "Validating Necklace source annotation against superTranscript coordinates"
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
log "Input preparation complete: $(wc -l < "$QUERY.fai") superTranscript sequences"
SORT_THREADS=1
MAP_THREADS=$((THREADS-3))
for species in D_simulans D_melanogaster; do
 log "========== $species: starting analysis =========="
 genome="$OUTROOT/tmp/$species.fa"
 if [[ $species == D_simulans ]]; then refgtf=$SIM_GTF; else refgtf=$MEL_GTF; fi
 log "$species | Validating reference annotation coordinates"
 python3 "$OUTROOT/tmp/validate.py" "$genome.fai" "$refgtf"
 bam="$OUTROOT/alignments/${species}_all.bam"
 log "$species | Aligning all superTranscripts with Minimap2 ($MAP_THREADS mapping threads); sorting BAM concurrently"
 minimap2 -ax splice:hq -u b --secondary=yes -N 5 --eqx -t "$MAP_THREADS" "$genome" "$QUERY" | samtools sort -@ "$SORT_THREADS" -m 1G -o "$bam"
 log "$species | Alignment and sorting complete; checking and indexing BAM"
 samtools quickcheck "$bam"; samtools index "$bam"
 log "$species | Computing mapping QC (flag counts, query coverage and identity)"
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

 log "$species | Extracting primary mapped alignments"
 primary="$OUTROOT/alignments/${species}_primary.bam"
 samtools view -b -F 2308 -o "$primary" "$bam"; samtools index "$primary"
 log "$species | Preparing JBrowse reference assembly and indexed reference annotation"
 track="$OUTROOT/jbrowse_tracks/$species"; mkdir -p "$track"
 cp "$genome" "$genome.fai" "$primary" "$primary.bai" "$track/"
 awk -F '\t' 'BEGIN{OFS="\t"} !/^#/ && NF==9' "$refgtf" | LC_ALL=C sort -k1,1 -k4,4n | bgzip -c > "$track/reference.gtf.gz"
 tabix -p gff "$track/reference.gtf.gz"
 log "$species | Generating BED12 and split alignment-block tracks"
 for mode in all primary; do
  b="$OUTROOT/alignments/${species}_${mode}.bam"
  bedtools bamtobed -bed12 -i "$b" > "$OUTROOT/bed/${species}_${mode}.bed12"
  bedtools bamtobed -split -i "$b" > "$OUTROOT/bed/${species}_${mode}_blocks.bed"
 done
 cp "$OUTROOT/bed/${species}_primary.bed12" "$track/"
 # Raw alignment-block depth, not expression; split at N and D.
 log "$species | Computing raw alignment depth (not RNA expression)"
 bedtools genomecov -ibam "$primary" -bg -split > "$OUTROOT/bigwig/${species}.bedgraph"
 if command -v bedGraphToBigWig >/dev/null; then
  cut -f1,2 "$genome.fai" > "$OUTROOT/tmp/$species.sizes"
  LC_ALL=C sort -k1,1 -k2,2n "$OUTROOT/bigwig/${species}.bedgraph" > "$OUTROOT/tmp/$species.sorted.bedgraph"
  if [[ -s $OUTROOT/tmp/$species.sorted.bedgraph ]]; then
   bedGraphToBigWig "$OUTROOT/tmp/$species.sorted.bedgraph" "$OUTROOT/tmp/$species.sizes" "$track/${species}_alignment_depth.bw"
  fi
 fi
 log "$species | Comparing primary alignment blocks with reference exons"
 # Candidate reference gene overlaps; keep raw rows for isoform-aware aggregation.
 awk -F '\t' '$3=="exon"' "$refgtf" > "$OUTROOT/tmp/$species.exons.gtf"
 bedtools intersect -a "$OUTROOT/bed/${species}_primary_blocks.bed" -b "$OUTROOT/tmp/$species.exons.gtf" -wo > "$OUTROOT/annotation/${species}_exon_overlaps.tsv"
 if [[ -n $DE_LIST ]]; then
  log "$species | Extracting requested DE-specific tracks"
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
   log "$species | Liftoff trial: distance=$d, coverage=0.5, identity=0.5, options=${opts[*]:-default}"
   log "$species | Liftoff detailed progress: $dest/liftoff.log"
   if (cd "$dest"; liftoff "$genome" "$QUERY" -g "$SOURCE" -o projected.annotation -u unmapped.txt -dir intermediate -p "$THREADS" -d "$d" -a 0.5 -s 0.5 "${opts[@]}") 2>&1 | tee "$dest/liftoff.log"; then
    log "$species | Liftoff distance=$d finished; validating projected annotation"
    # Liftoff preserves source format; detect from output rather than suffix.
    if grep -m1 -q '^##gff-version 3' "$dest/projected.annotation"; then ext=gff3; else ext=${SOURCE##*.}; fi
    cp "$dest/projected.annotation" "$dest/projected.$ext"
    if python3 "$OUTROOT/tmp/validate.py" "$genome.fai" "$dest/projected.$ext"; then
     cp "$dest/projected.$ext" "$track/Necklace_d${d}.$ext"
     if [[ $ext == gtf ]] && command -v gffread >/dev/null; then
      gffread -E -F --keep-genes "$dest/projected.$ext" -o "$track/Necklace_d${d}.converted.gff3" 2> "$dest/gffread.log" || { rm -f "$track/Necklace_d${d}.converted.gff3"; echo "GFF conversion failed; native projection retained"; }
     fi
     log "$species | Liftoff distance=$d: projection validated and exported"
     printf '%s\t%s\t%s\n' "$species" "$d" OK >> "$OUTROOT/qc/liftoff_status.tsv"
    else log "$species | Liftoff distance=$d: INVALID projection"; printf '%s\t%s\t%s\n' "$species" "$d" INVALID >> "$OUTROOT/qc/liftoff_status.tsv"; fi
   else log "$species | Liftoff distance=$d: FAILED; see $dest/liftoff.log"; printf '%s\t%s\t%s\n' "$species" "$d" FAILED >> "$OUTROOT/qc/liftoff_status.tsv"; fi
  done
 fi
 log "========== $species: analysis finished (check Liftoff trial status) =========="
done
printf 'Finished %s\nOutputs: %s\n' "$(date -Is)" "$OUTROOT"
if [[ $RUN_LIFTOFF == 1 ]] && grep -Eq 'FAILED|INVALID' "$OUTROOT/qc/liftoff_status.tsv"; then echo 'One or more Liftoff trials failed; inspect logs. Mapping outputs retained.'; exit 2; fi
