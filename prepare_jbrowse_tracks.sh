#!/usr/bin/env bash
# FAIR browser export, not a new mapping stage. Original PSL-to-GTF useful for
# localisation but lacks source block relationships. Retain native annotation,
# primary splice-aware BAM, direct graph block fragments and both audited
# Liftoff d2 projections separately. Graph blocks/fragments are not native exons
# and are distinct from the counted blocks used in limma. Source IDs preserved
# in direct track; Liftoff track uses reversible lookup and original_id attrs.
# Do not use old 165-gene Liftoff tracks. No inferred orthology.
# Output is a portable set of indexed files, not a configured JBrowse instance.
# Record SHA256, commands and versions; never modify existing inputs.
set -Eeuo pipefail
MAP=${MAP:-SuperTranscript_mapping_20261125_215819}
DIRECT=${DIRECT:-Necklace_CIGAR_projection}
SIM_RUN=${SIM_RUN:-Necklace_annotation_projection_20261125_224113}
MEL_RUN=${MEL_RUN:-Necklace_annotation_projection_20261125_225413}
OUT=${OUT:-Necklace_JBrowse_tracks}
CANDIDATES=${CANDIDATES:-}
for t in samtools bgzip tabix python3; do command -v "$t" >/dev/null; done
[[ ! -e "$OUT" ]] || { echo "Output exists: $OUT"; exit 1; }
mkdir -p "$OUT"
exec > >(tee -a "$OUT/export.log") 2>&1
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
{ date -Is; samtools --version; bgzip --version; tabix --version; printf 'MAP=%s\nDIRECT=%s\nSIM_RUN=%s\nMEL_RUN=%s\nCANDIDATES=%s\n' "$MAP" "$DIRECT" "$SIM_RUN" "$MEL_RUN" "$CANDIDATES"; sha256sum "$0"; } > "$OUT/provenance.txt"
pack_gff(){
 local src=$1 dst=$2
 # Preserve GFF3 header only; coordinate-sort records for tabix. Parent
 # relationships remain encoded in attributes even if file order changes.
 { printf '##gff-version 3\n'; awk -F '\t' '/^##FASTA/{exit} !/^#/ && NF==9' "$src" | LC_ALL=C sort -k1,1 -k4,4n -k5,5n; } | bgzip -c > "$dst"
 tabix -p gff "$dst"
 sha256sum "$src" >> "$OUT/provenance.txt"
}
for sp in D_simulans D_melanogaster; do
 run=$SIM_RUN
 [[ $sp != D_melanogaster ]] || run=$MEL_RUN
 log "$sp: exporting independent tracks"
 dest="$OUT/$sp"
 mkdir "$dest"
 ref="$MAP/jbrowse_tracks/$sp"
 for f in "$sp.fa" "$sp.fa.fai" "$sp"_primary.bam "$sp"_primary.bam.bai reference.gtf.gz reference.gtf.gz.tbi; do
  cp "$ref/$f" "$dest/"
  sha256sum "$ref/$f" >> "$OUT/provenance.txt"
 done
 for mode in whole blocks; do
  awk -F '\t' -v sp="$sp" -v mode="$mode" '$1==sp && $2==mode && $3==2 && $NF=="AUDITED"{ok=1} END{exit !ok}' "$run/projection_comparison.tsv"
  pack_gff "$run/liftoff/$sp/${mode}_d2/projected.gff3" "$dest/Necklace_Liftoff_${mode}_d2.gff3.gz"
 done
 cp "$run/prepared_inputs/superTranscript_ID_lookup.tsv" "$dest/"
 pack_gff "$DIRECT/$sp/Necklace_CIGAR_blocks.gff3" "$dest/Necklace_CIGAR_blocks.gff3.gz"
 if [[ -n "$CANDIDATES" ]]; then
  sha256sum "$CANDIDATES" >> "$OUT/provenance.txt"
  samtools view -b -N "$CANDIDATES" -o "$dest/Necklace_candidates.bam" "$dest/${sp}_primary.bam"
  samtools index "$dest/Necklace_candidates.bam"
 fi
 samtools quickcheck "$dest/${sp}_primary.bam"
done
log "Export complete: $OUT"
