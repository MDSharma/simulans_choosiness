#!/usr/bin/env bash
# Necklace -> genome annotation projection benchmark | M. D. Sharma | 2026-11
#
# SCIENTIFIC RATIONALE
# The original BLAT PSL -> GTF converter establishes genomic locations but
# labels alignment blocks as exons, loses Necklace graph block identities,
# and does not distinguish alternative loci sharing a transcript ID.
# Liftoff uses annotation hierarchy to project features; this is complementary
# to direct splice-aware Minimap2, not necessarily more accurate. A Necklace
# superTranscript is a composite intron-free sequence, not a natural isoform.
# SuperDuper.gtf has one whole-sequence exon for counting; SuperDuper.gff
# has graph-derived blocks, not necessarily biological exons.
# In the 2026-11-25 run Liftoff extracted just 165/15319 sequences (163 REL-)
# with both d=2 and d=50, whereas direct Minimap2 mapped >97% to both species.
# Test explicit, sanitised whole vs graph-block hierarchies at fixed thresholds
# and two distance factors. This addresses extraction before mapping quality.
# Do not assume a projected block is a true exon, or agreement proves orthology.
#
# FAIR PRINCIPLES
# Findable: labelled trial directories, timestamped progress, summary TSV.
# Accessible/interoperable: FASTA, GFF3, TSV, original ID lookup, open tools.
# Reusable: input/script SHA256, versions, command logs, exact parameters,
# retained intermediates and failure diagnostics. Original input files and
# completed Minimap2 output are never modified. Archive the lookup with tracks.
# FAIR-oriented design, not certification; repository access/licensing and
# data deposition remain the project owner's decisions.
# Citations: Davidson et al. 2017 doi:10.1186/s13059-017-1284-1;
# Davidson & Oshlack 2018 doi:10.1093/gigascience/giy045;
# Shumate & Salzberg 2021 doi:10.1093/bioinformatics/btaa1016;
# Li 2018 doi:10.1093/bioinformatics/bty191.
set -Eeuo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKDIR=$(realpath "${WORKDIR:-$PWD}")
cd "$WORKDIR"
THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
OUTROOT=$(realpath -m "${OUTROOT:-Necklace_annotation_projection_$(date +%Y%m%d_%H%M%S)}")
MAPPING_RUN=$(realpath "${MAPPING_RUN:-SuperTranscript_mapping_20261125_215819}")
SPECIES=${SPECIES:-D_simulans}
LIFT_DISTANCES=${LIFT_DISTANCES:-'2 50'}
FASTA=$(realpath "${SUPERTRANSCRIPTS:-superTranscriptome/SuperDuper.fasta}")
BLOCKS=$(realpath "${NECKLACE_BLOCKS:-superTranscriptome/SuperDuper.gff}")
[[ $THREADS =~ ^[1-9][0-9]*$ ]] || { echo "Invalid THREADS"; exit 1; }
[[ ! -e $OUTROOT ]] || { echo "Output exists: $OUTROOT; choose a new OUTROOT"; exit 1; }
mkdir -p "$OUTROOT"/{logs,qc,liftoff,targets}
exec > >(tee -a "$OUTROOT/logs/pipeline.log") 2>&1
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
trap 'rc=$?; log "FAILED line $LINENO (exit $rc): $BASH_COMMAND" >&2' ERR
log "Starting annotation projection benchmark; CPUs=$THREADS; species=$SPECIES"
log "Using previous mapping: $MAPPING_RUN; results: $OUTROOT"
for t in python3 liftoff samtools sha256sum; do
 command -v "$t" >/dev/null || { log "Missing dependency: $t"; exit 1; }
done
python3 -c 'import gffutils' || { log "Python environment needs gffutils"; exit 1; }
[[ -s $FASTA && -s $BLOCKS && -s $MAPPING_RUN/tmp/query.fa ]] || { log "Missing source or completed mapping inputs"; exit 1; }
# Require identical source sequence file to that used in direct mapping.
cmp -s "$FASTA" "$MAPPING_RUN/tmp/query.fa" || { log "Source FASTA differs from original mapping input; aborting"; exit 1; }
{
 date -Is
 printf 'THREADS=%s\nSPECIES=%s\nLIFT_DISTANCES=%s\nMAPPING_RUN=%s\n' "$THREADS" "$SPECIES" "$LIFT_DISTANCES" "$MAPPING_RUN"
 python3 --version
 liftoff -V
 samtools --version
 sha256sum "$FASTA" "$BLOCKS" "$SCRIPT_DIR/prepare_necklace_annotation.py" "$SCRIPT_DIR/evaluate_liftoff_projection.sh"
 if command -v git >/dev/null; then git -C "$SCRIPT_DIR" rev-parse HEAD 2>/dev/null || true; fi
} > "$OUTROOT/logs/provenance.txt"
log "Preparing reversible identifiers and validating both explicit GFF3 hierarchies"
python3 -u "$SCRIPT_DIR/prepare_necklace_annotation.py" --fasta "$FASTA" --blocks "$BLOCKS" --out "$OUTROOT/prepared_inputs"
prep="$OUTROOT/prepared_inputs"
query="$prep/SuperDuper.standardised.fasta"
samtools faidx "$query"
summary="$OUTROOT/projection_comparison.tsv"
printf 'species\tannotation\tdistance\texpected_genes\textracted_genes\tprojected_genes\tunmapped_list_entries\tpartial_flagged\tlow_identity_flagged\texpected_exon_features\tprojected_exon_features\tprimary_span_overlap\tprimary_span_comparable\tstatus\n' > "$summary"
# Each trial gets its own copy of source annotation and writable work directory.
# The genome copy and its default Minimap2 index may be reused within a species.
# No custom splice preset inside Liftoff in this controlled first experiment.
failures=0
for species in $SPECIES; do
 case "$species" in D_simulans|D_melanogaster) ;; *) log "Unknown species $species"; exit 1;; esac
 genome="$OUTROOT/targets/$species.fa"
 bed="$MAPPING_RUN/bed/${species}_primary.bed12"
 [[ -s $MAPPING_RUN/tmp/$species.fa && -s $bed ]] || { log "Missing prior mapping inputs for $species"; exit 1; }
 log "$species | Copying target FASTA for writable Liftoff indexes"
 cp "$MAPPING_RUN/tmp/$species.fa" "$genome"
 samtools faidx "$genome"
 sha256sum "$genome" "$bed" >> "$OUTROOT/logs/provenance.txt"
 for mode in whole blocks; do
  for d in $LIFT_DISTANCES; do
   [[ $d =~ ^[0-9]+([.][0-9]+)?$ ]] || { log "Invalid distance $d"; exit 1; }
   dest="$OUTROOT/liftoff/$species/${mode}_d$d"
   mkdir -p "$dest"
   cp "$prep/SuperDuper.$mode.gff3" "$dest/source.gff3"
   log "$species | $mode d=$d | Liftoff coverage=0.5 identity=0.5 (default internal alignment)"
   rc=0
   (cd "$dest"; set -x; liftoff "$genome" "$query" -g source.gff3 -o projected.gff3 -u unmapped.txt -dir intermediate -p "$THREADS" -d "$d" -a 0.5 -s 0.5) 2>&1 | tee "$dest/liftoff.log" || rc=$?
   log "$species | $mode d=$d | auditing extracted IDs, projected features and primary alignment loci"
   if python3 - "$prep/superTranscript_ID_lookup.tsv" "$dest" "$genome.fai" "$bed" "$summary" "$species" "$mode" "$d" "$rc" <<'PY'
import csv, pathlib, sys
from urllib.parse import unquote
lut, dest, fai, bed, summary, species, mode, distance, rc = sys.argv[1:]
dest = pathlib.Path(dest)
with open(lut) as f:
    lookup = {r["gene_id"]: r for r in csv.DictReader(f, delimiter="\t")}
expected = set(lookup)
expected_exons = sum(int(r["block_count"]) if mode == "blocks" else 1 for r in lookup.values())
extracted_path = dest / "intermediate/reference_all_genes.fa"
extracted = []
if extracted_path.exists():
    with extracted_path.open() as f:
        extracted = [l[1:].split()[0] for l in f if l.startswith(">")]
actual = set(extracted)
with (dest / "missing_extracted_gene_ids.tsv").open("w") as f:
    f.write("gene_id\toriginal_id\n")
    for gid in sorted(expected - actual):
        f.write(f"{gid}\t{lookup[gid]['original_id']}\n")
lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
primary = {}
with open(bed) as f:
    for l in f:
        a = l.rstrip("\n").split("\t")
        primary.setdefault(a[3], []).append((a[0], int(a[1]), int(a[2]), a[5]))
genes, exons, invalid = {}, set(), 0
projected = dest / "projected.gff3"
if projected.exists():
    with projected.open() as f:
        for l in f:
            if l.startswith("#") or not l.strip():
                continue
            a = l.rstrip("\n").split("\t")
            try:
                if len(a) != 9 or not 1 <= int(a[3]) <= int(a[4]) <= lengths[a[0]]:
                    raise ValueError()
                attr = dict(x.split("=", 1) for x in a[8].strip(";").split(";") if "=" in x)
                fid = attr["ID"]
                if a[2] == "gene":
                    if fid in genes or fid not in expected:
                        raise ValueError()
                    genes[fid] = (a[0], int(a[3])-1, int(a[4]), a[6], attr)
                elif a[2] == "exon":
                    exons.add(fid)
            except (ValueError, KeyError):
                invalid += 1
partial = low = comparable = overlap = 0
with (dest / "per_gene_projection.tsv").open("w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["original_id", "gene_id", "extracted", "projected", "contig", "start_1based", "end", "strand", "coverage", "sequence_ID", "partial_mapping", "low_identity", "primary_span_overlap"])
    for gid, r in lookup.items():
        name = r["original_id"]
        if gid not in genes:
            w.writerow([name, gid, int(gid in actual), 0] + ["NA"] * 9)
            continue
        contig, start, end, strand, a = genes[gid]
        pf = a.get("partial_mapping", "NA")
        lf = a.get("low_identity", "NA")
        partial += pf.lower() in ("true", "1", "yes")
        low += lf.lower() in ("true", "1", "yes")
        agreement = "NA"
        if name in primary:
            comparable += 1
            # A coarse locus comparison, NOT exon/block concordance.
            match = any(c == contig and max(s, start) < min(e, end)
                        for c, s, e, st in primary[name])
            overlap += match
            agreement = int(match)
        w.writerow([name, gid, int(gid in actual), 1, contig, start+1, end, strand,
                    a.get("coverage", "NA"), a.get("sequence_ID", "NA"), pf, lf, agreement])
unmapped = dest / "unmapped.txt"
unmapped_n = sum(bool(l.strip()) for l in unmapped.open()) if unmapped.exists() else "NA"
if rc != "0":
    status = "LIFTOFF_FAILED"
elif actual != expected or len(extracted) != len(expected):
    status = "EXTRACTION_MISMATCH"
elif invalid or not projected.exists() or not unmapped.exists():
    status = "INVALID_OUTPUT"
else:
    status = "AUDITED"  # not evidence of accurate biological projection
with open(summary, "a") as f:
    f.write("\t".join(map(str, [species, mode, distance, len(expected), len(extracted),
             len(genes), unmapped_n, partial, low, expected_exons, len(exons),
             overlap, comparable, status])) + "\n")
print(f"Extracted {len(extracted)}/{len(expected)}; projected {len(genes)} genes and {len(exons)} distinct exon IDs; primary span overlap {overlap}/{comparable}; {status}", flush=True)
if status != "AUDITED":
    sys.exit(2)
PY
   then
    log "$species | $mode d=$d | audit complete"
   else
    log "$species | $mode d=$d | audit failed; inspect trial outputs"
    failures=$((failures+1))
   fi
  done
 done
done
log "Benchmark complete. Summary: $summary"
cat "$summary"
if (( failures > 0 )); then log "$failures trial(s) failed audit; results retained"; exit 2; fi
log "All trials passed technical audit; inspect projection quality before using tracks."
