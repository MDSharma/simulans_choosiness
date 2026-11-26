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
OUT=${OUT:-Necklace_JBrowse_portable}
BLOCKS=${BLOCKS:-Tan_original/DE_expression_analysis/blocks.gtf}
MASTER=${MASTER:-Necklace_integrated_annotation_native/master_annotation.tsv}
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
 for f in "$sp.fa" "$sp.fa.fai" "$sp"_primary.bam "$sp"_primary.bam.bai reference.gtf.gz reference.gtf.gz.tbi "$sp"_primary.bed12; do
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
mkdir -p "$OUT/tables" "$OUT/source" "$OUT/scripts"
cp "$MASTER" "$OUT/tables/master_annotation.tsv"
cp "$MAP/tmp/query.fa" "$OUT/source/SuperDuper.fasta"
cp "$MEL_RUN/prepared_inputs/SuperDuper.blocks.gff3" "$OUT/source/"
cp "$MEL_RUN/prepared_inputs/superTranscript_ID_lookup.tsv" "$OUT/source/"
cp "$DIRECT/summary.tsv" "$OUT/tables/direct_projection_summary.tsv"
cp Necklace_projection_comparison/block_comparison_summary.tsv "$OUT/tables/"
cp Necklace_projection_comparison/liftoff_shared_blocks.tsv "$OUT/tables/"
for sp in D_simulans D_melanogaster; do
 mkdir -p "$OUT/tables/$sp"
 cp Necklace_projection_comparison/"$sp"/{gene_summary.tsv,reference_gene_candidates.tsv} "$OUT/tables/$sp/"
 cp "$MAP/qc/${sp}_alignments.tsv" "$OUT/tables/$sp/"
done
for f in prepare_jbrowse_tracks.sh project_necklace_blocks.py project_counted_blocks.py configure_portable_jbrowse.py integrate_superTranscript_annotations.py summarise_candidate_superTranscripts.py; do
 cp "$f" "$OUT/scripts/"
done
if [[ -f "$BLOCKS" ]]; then
 cp "$BLOCKS" "$OUT/source/blocks.gtf"
 python3 -u project_counted_blocks.py --mapping-run "$MAP" --blocks "$BLOCKS" --out "$OUT/tables/counted_blocks"
 for sp in D_simulans D_melanogaster; do
  pack_gff "$OUT/tables/counted_blocks/$sp.gff3" "$OUT/$sp/Necklace_counted_blocks.gff3.gz"
 done
else
 log "No counted block GTF at $BLOCKS; optional track omitted"
fi
python3 configure_portable_jbrowse.py "$OUT"
# Optional offline browser runtime: copy a previously extracted JBrowse 2 web
# release (its index.html and assets); do not silently download software.
if [[ -n "${JBROWSE_WEB_DIR:-}" ]]; then
 [[ -f "$JBROWSE_WEB_DIR/index.html" ]]
 mkdir "$OUT/web"
 cp -a "$JBROWSE_WEB_DIR/." "$OUT/web/"
 cp "$OUT/config.json" "$OUT/web/config.json"
 # Relative track URIs resolve against config location. Serve from root,
 # using web/index.html?config=../config.json, not web/config.json.
 rm "$OUT/web/config.json"
fi
(cd "$OUT"; find . -type f ! -name SHA256SUMS ! -name export.log -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > SHA256SUMS)
log "Export complete: $OUT"
