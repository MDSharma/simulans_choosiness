
#!/usr/bin/env bash
set -euo pipefail

#######################################################################
#  SUPERTRANSCRIPT → GENOME MAPPING PIPELINE
#  With Liftoff annotation + JBrowse 2 Desktop-ready outputs
#
#  --------------------------------------------------------------------
#  FAIR-COMPLIANT PIPELINE DESCRIPTION
#  --------------------------------------------------------------------
#  This pipeline follows the FAIR principles (Findable, Accessible,
#  Interoperable, Reusable) to ensure all generated data and metadata
#  can be reused, inspected, and understood long-term.
#
#  FINDABLE:
#    - Output files are organised into structured directories (alignments/,
#      liftoff/, bigwig/, bed/, jbrowse_tracks/, tmp/, logs/) with clear names.
#    - BAM/GFF3 headers and file names carry sufficient metadata to be
#      discovered by IGV/JBrowse and downstream tools.
#
#  ACCESSIBLE:
#    - Outputs use open, widely supported formats:
#        BAM/BAI (alignment), BED (intervals), BigWig (coverage),
#        GFF3 (annotations), FASTA (sequences), TSV (BLAST hits).
#    - Uses open-source, community-maintained tools (minimap2, samtools,
#      bedtools, Liftoff, deepTools).
#
#  INTEROPERABLE:
#    - Conforms to standard schemas/specs (SAM/BAM, BED, BigWig, GFF3).
#    - Liftoff maintains the same transcript/gene identifiers from Necklace,
#      easing cross-study integration and comparison.
#    - File/folder conventions interoperate with Nextflow/Snakemake,
#      R/Bioconductor, Python libraries, and common viewers.
#
#  REUSABLE:
#    - Deterministic, transparent steps with explicit tool versions logged.
#    - Self-contained, modular script with comments and clear defaults.
#    - The jbrowse_tracks/ folder can be archived/shared (e.g., Zenodo/OSF)
#      for long-term reuse and citation. Liftoff unmapped reports support QC.
#
#  --------------------------------------------------------------------
#  BENEFITS OF THIS APPROACH (vs PSL→GTF conversion alone)
#  --------------------------------------------------------------------
#  1. Minimap2 generates accurate, splice-aware alignments of long
#     superTranscripts, robust across species divergence.
#  2. BAM tracks preserve rich alignment metadata for IGV/JBrowse:
#     mismatches, CIGAR, splice junctions, mapping qualities.
#  3. Using the Necklace GTF with Liftoff produces biologically meaningful,
#     genome-specific annotation (GFF3) comparable across species.
#  4. JBrowse 2 Desktop-ready organisation: drag-and-drop tracks; no JSON
#     editing or CLI configuration required.
#
#######################################################################

#######################################################################
#  AUTHORSHIP, ORCID & CITATION METADATA
#
#  Pipeline Author(s):
#    Name: Sharma, M D
#    ORCID: https://orcid.org/0000-0002-9957-3153
#
#  Citation Instructions:
#    If you use this pipeline, please cite:
#      Sharma, M D. (2026). "SuperTranscript → Genome Mapping and
#      Annotation Pipeline." Version 1.0. DOI: [Zenodo DOI will be made available later]
#      Repository: https://github.com/MDSharma/simulans_choosiness
#
#  Software Citations:
#    Li H. (2018) Minimap2. Bioinformatics 34:3094–3100.
#    Danecek P. (2021) SAMtools/BCFtools. GigaScience 10:giab008.
#    Quinlan AR. (2010) BEDTools. Bioinformatics 26:841–842.
#    Shumate A., Salzberg SL. (2021) Liftoff. Bioinformatics 37:1639–1643.
#    Ramírez F. (2016) deepTools2. NAR 44:W160–W165.
#    Buels R. (2016) JBrowse. Genome Biology 17:66.
#    Necklace / SuperTranscripts:
#      Davidson NM. (2017) Genome Research 27:2018–2029.
#
#  License:
#    Released under the [MIT / GPL-3.0 / CC-BY 4.0] license.
#
#######################################################################

#######################################################################
#  DEPENDENCY CHECK
#  Ensures required tools are available; hints for installation provided.
#######################################################################
echo "[DEPCHK] Checking dependencies..."

REQUIRED_TOOLS=( samtools minimap2 bedtools bamCoverage liftoff awk )
OPTIONAL_TOOLS=( makeblastdb blastn )

MISSING=0

check_tool() {
  local t="$1"
  if ! command -v "$t" >/dev/null 2>&1; then
    echo "❌  Missing: $t"
    case "$t" in
      samtools)
        echo "    Install Samtools: https://www.htslib.org/download/ | mamba install -c bioconda samtools" ;;
      minimap2)
        echo "    Install Minimap2: https://github.com/lh3/minimap2 | mamba install -c bioconda minimap2" ;;
      bedtools)
        echo "    Install BEDTools: https://github.com/arq5x/bedtools2 | mamba install -c bioconda bedtools" ;;
      bamCoverage)
        echo "    Install deepTools (bamCoverage): https://deeptools.readthedocs.io | mamba install -c bioconda deeptools" ;;
      liftoff)
        echo "    Install Liftoff: https://github.com/agshumate/liftoff | mamba install -c bioconda liftoff" ;;
      makeblastdb|blastn)
        echo "    Install BLAST+: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ | mamba install -c bioconda blast" ;;
      awk)
        echo "    Install awk (gawk): via system package manager (e.g., apt, yum, brew)." ;;
    esac
    MISSING=1
  else
    echo "✔  Found: $t ($(command -v "$t"))"
  fi
}

echo "[DEPCHK] Required:"
for t in "${REQUIRED_TOOLS[@]}"; do check_tool "$t"; done

echo "[DEPCHK] Optional (for BLAST step):"
for t in "${OPTIONAL_TOOLS[@]}"; do
  if command -v "$t" >/dev/null 2>&1; then
    echo "✔  Found: $t"
  else
    echo "⚠  Optional missing: $t (BLAST annotation will be skipped if DBs not provided)"
  fi
done

if [[ $MISSING -eq 1 ]]; then
  echo "❌ One or more required tools are missing. Install them and re-run."
  exit 1
fi

echo "------------------------------------------------------------"

#######################################################################
#  USER-DEFINED VARIABLES
#######################################################################

# SuperTranscriptome FASTA and Necklace GTF
SUPERTRANSCRIPTS="superTranscriptome.fasta"
SUPERTRANSCRIPT_GTF="superTranscriptome.gtf"   # <- Necklace-produced GTF/GFF3

# Differentially expressed superTranscript IDs (one per line; must match FASTA headers)
DE_LIST="DE_transcript_ids.txt"

# Reference genomes (target species)
GENOME_A="speciesA.fasta"
GENOME_B="speciesB.fasta"

# Optional: CDS/Transcript FASTAs to run BLAST annotations against (if present, BLAST will run)
GENESET_A="speciesA_cds.fasta"
GENESET_B="speciesB_cds.fasta"

# Threads
THREADS=8

#######################################################################
#  INPUT SANITY CHECKS
#######################################################################
echo "[INPUT] Verifying input files exist..."

for f in "$SUPERTRANSCRIPTS" "$SUPERTRANSCRIPT_GTF" "$DE_LIST" "$GENOME_A" "$GENOME_B"; do
  if [[ ! -f "$f" ]]; then
    echo "❌ Missing input: $f"
    exit 1
  fi
done
echo "✔ All required inputs present."

#######################################################################
#  OUTPUT STRUCTURE
#######################################################################
RUNSTAMP=$(date +"%Y%m%d_%H%M%S")
OUTROOT="run_${RUNSTAMP}"
mkdir -p "$OUTROOT"/{alignments,bed,bigwig,liftoff,jbrowse_tracks,blast,tmp,logs}

echo "[VERSIONS] Capturing tool versions..." | tee "logs/versions.txt"
{
  echo "Date: $(date -Is)"
  echo -n "samtools "; samtools --version | head -n1
  echo -n "minimap2 "; minimap2 --version
  echo -n "bedtools "; bedtools --version 2>/dev/null || true
  echo -n "bamCoverage "; bamCoverage --version 2>/dev/null || true
  echo -n "liftoff "; liftoff --version 2>/dev/null || true
  echo -n "blastn "; (blastn -version 2>/dev/null | head -n1) || true
} >> "logs/versions.txt"

#######################################################################
#  1) EXTRACT DE superTranscripts
#######################################################################
echo "[1/8] Extracting DE superTranscripts..." | tee -a "logs/steps.log"

samtools faidx "$SUPERTRANSCRIPTS"
xargs samtools faidx "$SUPERTRANSCRIPTS" < "$DE_LIST" > "$OUTROOT/tmp/DE_superTranscripts.fasta"

# Quick count
DE_N=$(grep -c '^>' "$OUTROOT/tmp/DE_superTranscripts.fasta" || true)
echo "   Extracted $DE_N sequences." | tee -a "logs/steps.log"

#######################################################################
#  2) ALIGN TO BOTH GENOMES (minimap2 → BAM)
#######################################################################
echo "[2/8] Aligning DE superTranscripts to Genome A..." | tee -a "logs/steps.log"
minimap2 -ax splice -uf --secondary=no -t "$THREADS" "$GENOME_A" "$OUTROOT/tmp/DE_superTranscripts.fasta" \
  | samtools sort -@ "$THREADS" -o "$OUTROOT/alignments/speciesA_DE.bam"
samtools index "$OUTROOT/alignments/speciesA_DE.bam"

echo "[2/8] Aligning DE superTranscripts to Genome B..." | tee -a "logs/steps.log"
minimap2 -ax splice -uf --secondary=no -t "$THREADS" "$GENOME_B" "$OUTROOT/tmp/DE_superTranscripts.fasta" \
  | samtools sort -@ "$THREADS" -o "$OUTROOT/alignments/speciesB_DE.bam"
samtools index "$OUTROOT/alignments/speciesB_DE.bam"

#######################################################################
#  3) BAM → BED (exons) and merged spans
#######################################################################
echo "[3/8] Converting BAM to BED + merged spans..." | tee -a "logs/steps.log"

bedtools bamtobed -i "$OUTROOT/alignments/speciesA_DE.bam" > "$OUTROOT/bed/speciesA_DE.bed"
bedtools bamtobed -i "$OUTROOT/alignments/speciesB_DE.bam" > "$OUTROOT/bed/speciesB_DE.bed"

bedtools sort -i "$OUTROOT/bed/speciesA_DE.bed" | bedtools merge > "$OUTROOT/bed/speciesA_DE_merged.bed"
bedtools sort -i "$OUTROOT/bed/speciesB_DE.bed" | bedtools merge > "$OUTROOT/bed/speciesB_DE_merged.bed"

#######################################################################
#  4) COVERAGE TRACKS (BigWig via deepTools)
#######################################################################
echo "[4/8] Building BigWig coverage (CPM-normalised)..." | tee -a "logs/steps.log"

bamCoverage -b "$OUTROOT/alignments/speciesA_DE.bam" -o "$OUTROOT/bigwig/speciesA_DE.bw" \
  --normalizeUsing CPM --binSize 1 --numberOfProcessors "$THREADS" \
  --ignoreDuplicates --minMappingQuality 1

bamCoverage -b "$OUTROOT/alignments/speciesB_DE.bam" -o "$OUTROOT/bigwig/speciesB_DE.bw" \
  --normalizeUsing CPM --binSize 1 --numberOfProcessors "$THREADS" \
  --ignoreDuplicates --minMappingQuality 1

#######################################################################
#  5) ANNOTATION LIFTOVER WITH LIFTOFF (using Necklace GTF)
#     NOTE: Liftoff syntax is: liftoff <target.fa> <reference.fa> -g <reference.gtf>
#######################################################################
echo "[5/8] Running Liftoff to lift Necklace annotation onto Genome A & B..." | tee -a "logs/steps.log"

liftoff "$GENOME_A" "$SUPERTRANSCRIPTS" -g "$SUPERTRANSCRIPT_GTF" \
  -o "$OUTROOT/liftoff/speciesA_liftoff.gff3" \
  -u "$OUTROOT/liftoff/speciesA_unmapped.txt" \
  -p "$THREADS"

liftoff "$GENOME_B" "$SUPERTRANSCRIPTS" -g "$SUPERTRANSCRIPT_GTF" \
  -o "$OUTROOT/liftoff/speciesB_liftoff.gff3" \
  -u "$OUTROOT/liftoff/speciesB_unmapped.txt" \
  -p "$THREADS"


#######################################################################
#  6) OPTIONAL BLAST ANNOTATION (runs only if reference CDS FASTAs exist)
#######################################################################
BLAST_DID_SOMETHING=0
if command -v makeblastdb >/dev/null 2>&1 && command -v blastn >/dev/null 2>&1; then
  echo "[6/8] BLAST+: Checking for provided gene sets..." | tee -a "logs/steps.log"
  if [[ -f "$GENESET_A" ]]; then
    echo "   Building BLAST DB for species A and running blastn..." | tee -a "logs/steps.log"
    makeblastdb -in "$GENESET_A" -dbtype nucl -out "$OUTROOT/blast/speciesA_cds_db"
    blastn -query "$OUTROOT/tmp/DE_superTranscripts.fasta" \
           -db "$OUTROOT/blast/speciesA_cds_db" \
           -out "$OUTROOT/blast/speciesA_blast.tsv" -outfmt 6 -num_threads "$THREADS"
    BLAST_DID_SOMETHING=1
  else
    echo "   ⚠ No GENESET_A ($GENESET_A); skipping BLAST for A." | tee -a "logs/steps.log"
  fi
  if [[ -f "$GENESET_B" ]]; then
    echo "   Building BLAST DB for species B and running blastn..." | tee -a "logs/steps.log"
    makeblastdb -in "$GENESET_B" -dbtype nucl -out "$OUTROOT/blast/speciesB_cds_db"
    blastn -query "$OUTROOT/tmp/DE_superTranscripts.fasta" \
           -db "$OUTROOT/blast/speciesB_cds_db" \
           -out "$OUTROOT/blast/speciesB_blast.tsv" -outfmt 6 -num_threads "$THREADS"
    BLAST_DID_SOMETHING=1
  else
    echo "   ⚠ No GENESET_B ($GENESET_B); skipping BLAST for B." | tee -a "logs/steps.log"
  fi
else
  echo "[6/8] ⚠ BLAST+ not installed; skipping BLAST annotation." | tee -a "logs/steps.log"
fi

#######################################################################
#  7) PREPARE JBROWSE 2 DESKTOP EXPORT FOLDER
#######################################################################
echo "[7/8] Preparing jbrowse_tracks/ for drag-and-drop..." | tee -a "logs/steps.log"

# Core alignment + index
cp "$OUTROOT/alignments/speciesA_DE.bam" "$OUTROOT/jbrowse_tracks/"
cp "$OUTROOT/alignments/speciesA_DE.bam.bai" "$OUTROOT/jbrowse_tracks/"
cp "$OUTROOT/alignments/speciesB_DE.bam" "$OUTROOT/jbrowse_tracks/"
cp "$OUTROOT/alignments/speciesB_DE.bam.bai" "$OUTROOT/jbrowse_tracks/"

# Coverage
cp "$OUTROOT/bigwig/speciesA_DE.bw" "$OUTROOT/jbrowse_tracks/"
cp "$OUTROOT/bigwig/speciesB_DE.bw" "$OUTROOT/jbrowse_tracks/"

# BEDs (fine-grained + merged spans)
cp "$OUTROOT/bed/"speciesA_DE*.bed "$OUTROOT/jbrowse_tracks/" || true
cp "$OUTROOT/bed/"speciesB_DE*.bed "$OUTROOT/jbrowse_tracks/" || true

# Liftoff annotations
cp "$OUTROOT/liftoff/speciesA_liftoff.gff3" "$OUTROOT/jbrowse_tracks/"
cp "$OUTROOT/liftoff/speciesB_liftoff.gff3" "$OUTROOT/jbrowse_tracks/"

# (Optional) BLAST results (not a track, but packaged for convenience)
if [[ $BLAST_DID_SOMETHING -eq 1 ]]; then
  cp "$OUTROOT/blast/"species*_blast.tsv "$OUTROOT/jbrowse_tracks/" || true
fi

# (Optional) Also provide assemblies for convenience (users can load them in JBrowse Desktop)
cp "$GENOME_A" "$OUTROOT/jbrowse_tracks/" || true
cp "$GENOME_B" "$OUTROOT/jbrowse_tracks/" || true

#######################################################################
#  8) SUMMARY & NEXT STEPS
#######################################################################
echo
echo "==============================================================="
echo " PIPELINE COMPLETE — OUTPUT ROOT: $OUTROOT"
echo "---------------------------------------------------------------"
echo "  JBrowse 2 Desktop: drag the following into the viewer:"
echo "   - jbrowse_tracks/speciesA_DE.bam (+ .bai)"
echo "   - jbrowse_tracks/speciesB_DE.bam (+ .bai)"
echo "   - jbrowse_tracks/speciesA_DE.bw"
echo "   - jbrowse_tracks/speciesB_DE.bw"
echo "   - jbrowse_tracks/speciesA_liftoff.gff3"
echo "   - jbrowse_tracks/speciesB_liftoff.gff3"
echo "   - jbrowse_tracks/speciesA_DE.bed / speciesA_DE_merged.bed"
echo "   - jbrowse_tracks/speciesB_DE.bed / speciesB_DE_merged.bed"
echo "  (Assemblies copied for convenience: ${GENOME_A##*/}, ${GENOME_B##*/})"
echo "---------------------------------------------------------------"
echo " Logs:           logs/ (versions.txt, steps.log)"
echo " Alignments:     $OUTROOT/alignments/"
echo " Coverage:       $OUTROOT/bigwig/"
echo " Intervals:      $OUTROOT/bed/"
echo " Liftoff ann.:   $OUTROOT/liftoff/"
echo " BLAST:          $OUTROOT/blast/ (if provided)"
echo " Temp FASTAs:    $OUTROOT/tmp/"
echo "==============================================================="
echo
