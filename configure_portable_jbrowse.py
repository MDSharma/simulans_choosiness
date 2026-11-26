#!/usr/bin/env python3
"""Generate a portable JBrowse 2 configuration with relative paths only."""
import json,sys
from pathlib import Path
root=Path(sys.argv[1])
assemblies=[]; tracks=[]
for sp,label in [("D_simulans","D. simulans"),("D_melanogaster","D. melanogaster")]:
    def loc(name):
        path=f"{sp}/{name}"
        if not (root/path).is_file(): raise FileNotFoundError(path)
        return {"uri":path,"locationType":"UriLocation"}
    assemblies.append({"name":sp,"displayName":label,"sequence":{"type":"ReferenceSequenceTrack","trackId":sp+"_reference","adapter":{"type":"IndexedFastaAdapter","fastaLocation":loc(sp+".fa"),"faiLocation":loc(sp+".fa.fai")}}})
    tracks.append({"type":"AlignmentsTrack","trackId":sp+"_primary","name":"Minimap2 superTranscript alignments","assemblyNames":[sp],"adapter":{"type":"BamAdapter","bamLocation":loc(sp+"_primary.bam"),"index":{"location":loc(sp+"_primary.bam.bai"),"indexType":"BAI"}}})
    for name,title,kind in [("reference.gtf.gz","Native reference genes","Gtf"),("Necklace_Liftoff_whole_d2.gff3.gz","Liftoff whole superTranscripts","Gff3"),("Necklace_Liftoff_blocks_d2.gff3.gz","Liftoff graph blocks","Gff3"),("Necklace_CIGAR_blocks.gff3.gz","CIGAR graph block fragments","Gff3"),("Necklace_counted_blocks.gff3.gz","Recovered blocks.gtf (structural only)","Gff3")]:
        if name.startswith("Necklace_counted") and not (root/sp/name).exists(): continue
        tracks.append({"type":"FeatureTrack","trackId":sp+"_"+name.split(".")[0],"name":title,"assemblyNames":[sp],"adapter":{"type":kind+"TabixAdapter",kind.lower()+"GzLocation":loc(name),"index":{"location":loc(name+".tbi"),"indexType":"TBI"}}})
(root/"config.json").write_text(json.dumps({"assemblies":assemblies,"tracks":tracks},indent=2)+"\n")
(root/"README.md").write_text("""# Necklace superTranscript genome browser handover

Two assemblies, independent structural projections, all native gene candidates.
All paths in config.json are relative to this folder. No Athena paths needed.

## Open the browser

Use JBrowse 2 Desktop or an existing JBrowse 2 web installation. For a web
installation copied into web/ at packaging, serve this entire folder over HTTP
and open web/index.html?config=../config.json. Do not open through file://.
Use a server supporting HTTP byte ranges (e.g. Node http-server); from this
folder, run npx http-server . -p 8080 -c-1 and visit
http://localhost:8080/web/index.html?config=../config.json .
If web/ is absent, the folder contains all data and configuration but not the
browser application. Install JBrowse 2 separately or copy its web release into
web/. No external sequence or track data is needed.

Select either assembly in a linear genome view and open the desired tracks.
BED12 is supplied for optional manual import. Search native genes using the
annotation tables; no text search index is built.

## Interpretation and handover

* Native annotation is the reference GTF. Minimap2 BAM shows primary alignments.
* Liftoff whole and graph-block tracks use audited d=2 projections. IDs are
  linked to original Necklace IDs in the lookup table.
* CIGAR graph-block fragments are not biological exons; alignment orientation
  need not equal transcriptional strand.
* Recovered blocks.gtf projections are structural only. Their source row IDs
  are not limma exon IDs. Confirm count matrix, feature ordering and R
  make.unique() identifiers before attaching statistical results.
* tables/master_annotation.tsv contains all 15,319 source IDs; gene candidates
  are overlaps, not proof of orthology. Preserve ambiguous assignments.
* Historical count matrices contain 142,593 and 121,757 rows respectively.
  The manuscript reports the latter; exact correspondence remains unverified.
* Check the final group coding (commented group definitions and plotting labels
  in Limma_alt_splice.R disagree).
* Integrate the two Trinotate reports separately using the supplied script;
  the original pipeline inputs and full scripts are maintained in the GitHub
  repository MDSharma/simulans_choosiness. Statistical and functional annotation
  integration is intentionally left to the analyst. The supplied scripts are
  copies for provenance and are not a standalone rerunnable mapping pipeline.

Validate copied files with sha256sum -c SHA256SUMS from the package root.
""")
print("Wrote portable config.json and README.md",flush=True)
