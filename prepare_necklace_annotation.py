#!/usr/bin/env python3
"""Prepare explicit, reversible Necklace feature hierarchies for Liftoff.

M. D. Sharma | simulans_choosiness | 2026-11
Scientific rationale and FAIR provenance
----------------------------------------
Necklace/Lace produces a linear union of transcript sequence per inferred
gene cluster, not a conventional genomic transcript. SuperDuper.gtf describes
whole sequences for counting; SuperDuper.gff describes graph-derived blocks.
Neither graph blocks nor a composite superTranscript imply genuine exons or
a single natural isoform. Splice-junction-derived counts/blocks.gtf is a
different annotation and is not substituted here.

BLAT PSL -> GTF provides genomic alignment blocks, but the original simple
converter discarded source block identities and did not distinguish multiple
loci. Liftoff is annotation-aware; its suitability for intron-free source
sequences remains experimental. Direct splice-aware Minimap2 is the baseline.
On 2026-11-25, 15,319 queries gave 15,064 / 14,940 primary mappings to
D. simulans / D. melanogaster. Liftoff extracted only 165 source genes
(163 REL-), irrespective of distance factor, rather than testing all queries.
The precise parsing/inference defect is not established. Explicit hierarchy,
unique safe IDs and forward source orientation remove these confounders
together; this experiment does not isolate their individual effects.

FAIR: deterministic IDs in FASTA order, reversible TSV lookup, standard
FASTA/GFF3, SHA256 input and script fingerprints, software versions,
validation report, original inputs untouched. The derived gene/mRNA/exon
types are compatibility labels, not validated biological gene models.
2026-11-25: source transcript_id differs from seqid at line 139020 in
SuperDuper.gff. Require gene_id/seqid agreement and valid coordinates,
but audit rather than reject transcript label differences; retain original
attributes in source_transcript_id_mismatches.tsv. No blocks are discarded.
References:
 Davidson et al. 2017 doi:10.1186/s13059-017-1284-1
 Davidson & Oshlack 2018 doi:10.1093/gigascience/giy045
 Shumate & Salzberg 2021 doi:10.1093/bioinformatics/btaa1016
 https://github.com/Oshlack/necklace/wiki/Output
 https://github.com/Oshlack/Lace/wiki/Usage-Documentation
"""
import argparse
import collections
import csv
import hashlib
import pathlib
import re
import sys
from urllib.parse import quote

def log(message):
    print(message, flush=True)

def fasta(path):
    name, chunks = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError("Empty FASTA header")
                name, chunks = header[0], []
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def prepare(args):
    import gffutils
    sequences = {}
    log("Reading original superTranscript FASTA")
    for name, seq in fasta(args.fasta):
        if name in sequences or not seq or not re.fullmatch(r"[A-Za-z]+", seq):
            raise ValueError(f"Duplicate, empty or malformed FASTA sequence: {name}")
        sequences[name] = seq
    if not sequences:
        raise ValueError("Empty FASTA")
    log("Reading original Lace graph blocks (GTF-style attributes in .gff)")
    blocks = collections.defaultdict(list)
    transcript_id_mismatches = []
    with args.blocks.open() as handle:
        for n, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            if len(f) != 9 or f[2] != "exon" or f[0] not in sequences:
                raise ValueError(f"Unexpected graph block at line {n}: {line[:120]}")
            start, end = int(f[3]), int(f[4])
            if not 1 <= start <= end <= len(sequences[f[0]]):
                raise ValueError(f"Out-of-bounds block at line {n}")
            for key in ("gene_id", "transcript_id"):
                m = re.search(r'(?:^|;\s*)' + key + r'\s+"([^"]+)"', f[8])
                if not m:
                    raise ValueError(f"Missing or malformed {key} at line {n}")
                if m.group(1) != f[0]:
                    if key == "gene_id":
                        raise ValueError(f"gene_id differs from sequence ID at line {n}: {m.group(1)!r} vs {f[0]!r}")
                    # Graph blocks are located by seqid and coordinates, not by
                    # legacy transcript labels. Preserve discrepancies for audit;
                    # the derived hierarchy uses new, unique transcript IDs.
                    transcript_id_mismatches.append((n, f[0], start, end, m.group(1), f[8]))
            blocks[f[0]].append((start, end))
    missing = set(sequences) - set(blocks)
    if missing:
        raise ValueError(f"{len(missing)} FASTA sequences have no blocks; examples: {sorted(missing)[:5]}")
    # Lace graph blocks are expected to partition the linear sequence.
    # Do not silently repair gaps, duplicate blocks or overlaps.
    for name, intervals in blocks.items():
        intervals.sort()
        next_start = 1
        for start, end in intervals:
            if start != next_start:
                raise ValueError(f"Gap/overlap in {name}: expected {next_start}, got {start}")
            next_start = end + 1
        if next_start != len(sequences[name]) + 1:
            raise ValueError(f"Blocks do not cover end of {name}")
    nblocks = sum(map(len, blocks.values()))
    log(f"Validated {len(sequences)} sequences and {nblocks} graph blocks")
    log(f"Source transcript_id mismatches: {len(transcript_id_mismatches)} rows; retaining original labels in audit TSV")
    args.out.mkdir(parents=True, exist_ok=False)
    with (args.out / "source_transcript_id_mismatches.tsv").open("w") as audit:
        writer = csv.writer(audit, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_line", "original_id", "start", "end", "source_transcript_id", "source_attributes"])
        writer.writerows(transcript_id_mismatches)
    with (args.out / "SuperDuper.standardised.fasta").open("w") as fa, \
         (args.out / "SuperDuper.whole.gff3").open("w") as whole, \
         (args.out / "SuperDuper.blocks.gff3").open("w") as block, \
         (args.out / "superTranscript_ID_lookup.tsv").open("w") as lut:
        writer = csv.writer(lut, delimiter="\t", lineterminator="\n")
        writer.writerow(["projection_id", "original_id", "gene_id", "transcript_id", "length", "block_count"])
        for out in (whole, block):
            out.write("##gff-version 3\n")
        for i, (name, seq) in enumerate(sequences.items(), 1):
            sid = f"ST{i:06d}"
            gid, tid = f"gene_{sid}", f"tx_{sid}"
            original = quote(name, safe=":._-")
            writer.writerow([sid, name, gid, tid, len(seq), len(blocks[name])])
            fa.write(f">{sid}\n")
            for j in range(0, len(seq), 80):
                fa.write(seq[j:j+80] + "\n")
            def emit(out, feature, start, end, attrs):
                out.write(f"{sid}\tNecklace\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}\n")
            for out, intervals, label in ((whole, [(1, len(seq))], "whole"), (block, blocks[name], "graph_block")):
                emit(out, "gene", 1, len(seq), f"ID={gid};Name={original};original_id={original}")
                emit(out, "mRNA", 1, len(seq), f"ID={tid};Parent={gid};original_id={original}")
                for j, (start, end) in enumerate(intervals, 1):
                    emit(out, "exon", start, end,
                         f"ID={label}_{sid}_{j:06d};Parent={tid};original_id={original};representation={label};block_number={j}")
    for label, expected in (("whole", len(sequences)), ("blocks", nblocks)):
        path = args.out / f"SuperDuper.{label}.gff3"
        log(f"Validating {label} hierarchy with gffutils {gffutils.__version__}")
        db = gffutils.create_db(str(path), dbfn=":memory:", merge_strategy="error",
                                disable_infer_genes=True, disable_infer_transcripts=True)
        genes = list(db.features_of_type("gene"))
        if len(genes) != len(sequences) or db.count_features_of_type("mRNA") != len(sequences) or db.count_features_of_type("exon") != expected:
            raise ValueError(f"Unexpected feature counts after gffutils parsing: {label}")
        for gene in genes:
            tx = list(db.children(gene, level=1, featuretype="mRNA"))
            if len(tx) != 1:
                raise ValueError(f"Invalid gene hierarchy: {gene.id}")
            exons = list(db.children(tx[0], level=1, featuretype="exon"))
            if not exons or any(e.seqid != gene.seqid for e in exons):
                raise ValueError(f"Invalid transcript hierarchy: {tx[0].id}")
        db.conn.close()
    with (args.out / "validation_summary.tsv").open("w") as out:
        out.write("metric\tvalue\n")
        for key, val in (("sequences", len(sequences)), ("graph_blocks", nblocks), ("whole_exons", len(sequences)),
                         ("source_transcript_id_mismatch_rows", len(transcript_id_mismatches)),
                         ("source_transcript_id_mismatch_sequences", len({r[1] for r in transcript_id_mismatches})),
                         ("python", sys.version.split()[0]), ("gffutils", gffutils.__version__),
                         ("fasta_path", args.fasta.resolve()), ("fasta_sha256", digest(args.fasta)),
                         ("blocks_path", args.blocks.resolve()), ("blocks_sha256", digest(args.blocks)),
                         ("script_sha256", digest(pathlib.Path(__file__))), ("status", "PASS")):
            out.write(f"{key}\t{val}\n")
    log(f"Preparation complete: {args.out}")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", type=pathlib.Path, default=pathlib.Path("superTranscriptome/SuperDuper.fasta"))
    p.add_argument("--blocks", type=pathlib.Path, default=pathlib.Path("superTranscriptome/SuperDuper.gff"))
    p.add_argument("--out", type=pathlib.Path, required=True)
    prepare(p.parse_args())

if __name__ == "__main__":
    main()
