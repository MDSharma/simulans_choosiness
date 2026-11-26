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
        tracks.append({"type":"FeatureTrack","trackId":sp+"_"+name.split(".")[0],"name":title,"assemblyNames":[sp],"adapter":{"type":kind+"TabixAdapter",("gffGzLocation" if kind=="Gff3" else "gtfGzLocation"):loc(name),"index":{"location":loc(name+".tbi"),"indexType":"TBI"}}})
config={"assemblies":assemblies,"tracks":tracks}
(root/"config.json").write_text(json.dumps(config,indent=2)+"\n")
# v4 Desktop must use explicit local paths on the recipient's machine.
# PowerShell requires no Python on Windows and resolves the package root
# at runtime rather than embedding Athena paths in the shared session.
(root/"Create_Windows_JBrowse.ps1").write_text(r'''$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
$config = Get-Content -LiteralPath (Join-Path $root 'config.json') -Raw | ConvertFrom-Json
function Convert-Location($node) {
    if ($null -eq $node) { return }
    if ($node -is [System.Array]) {
        foreach ($item in $node) { Convert-Location $item }
    } elseif ($node -is [System.Management.Automation.PSCustomObject]) {
        if ($node.locationType -eq 'UriLocation') {
            $path = Join-Path $root $node.uri
            if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Missing file: $path" }
            $path = (Resolve-Path -LiteralPath $path).Path
            $node.PSObject.Properties.Remove('uri')
            $node.locationType = 'LocalPathLocation'
            $node | Add-Member -NotePropertyName localPath -NotePropertyValue $path -Force
        } else {
            foreach ($prop in $node.PSObject.Properties) { Convert-Location $prop.Value }
        }
    }
}
Convert-Location $config
$out = Join-Path $root 'Necklace_Windows.jbrowse'
$json = $config | ConvertTo-Json -Depth 100
[System.IO.File]::WriteAllText($out, $json, (New-Object System.Text.UTF8Encoding($false)))
Write-Host "Created $out"
Write-Host "Open this file in JBrowse 2 Desktop. Regenerate after moving the folder."
''',encoding="utf-8")

(root/"README.md").write_text("# SuperTranscript mapping and annotation – female choosiness\n\nI've put together the genome browser files for both D. simulans and\nD. melanogaster using the final Necklace assembly (SuperDuper.fasta).\nThe idea is to be able to look at our superTranscripts alongside the\nexisting gene models in both species.\n\n## JBrowse\n\nI've used JBrowse 2 Desktop (v4.3.0) on Windows. Download the whole folder\nfrom OneDrive and keep the directory structure intact. Run the following\nfrom PowerShell in the downloaded folder:\n\n    powershell -NoProfile -ExecutionPolicy Bypass -File .\\Create_Windows_JBrowse.ps1\n\nOpen the resulting Necklace_Windows.jbrowse file in JBrowse. The script\nsets the local paths to the genome and track files, so run it again if\nyou move the folder. The config.json file uses relative paths and is\nintended for JBrowse Web.\n\nIf using the web version, with the JBrowse application installed in\nweb/, serve the package directory over HTTP with byte-range support\n(e.g. npx http-server . -p 8080 -c-1) and open\nhttp://localhost:8080/web/index.html?config=../config.json.\n\n## Tracks\n\nThere is a separate directory for each species, containing the reference\nFASTA and annotation, the primary Minimap2 alignments (BAM and BED12),\nand the Liftoff projections. I've retained the whole-superTranscript\nand graph-block Liftoff results separately (both using d=2).\n\nI've also projected the original Necklace graph blocks directly through\nthe Minimap2 alignments. This recovered more of the original blocks\nthan Liftoff, although the two methods generally agreed where both\nrecovered the same block. The fragments in this track should not be\ntreated as individual genomic exons.\n\nWhere available, I've included a separate projection of blocks.gtf\nfrom the original DE_expression_analysis directory. This is useful\nfor looking at the structure of the superTranscripts in the context\nof the differential exon usage analysis.\n\nThe BED12 files are available for manual import if needed. I haven't\nbuilt a text search index; the tables below can be used to find loci.\n\n## Annotation tables\n\nThe main file is tables/master_annotation.tsv, with one row for each\nof the 15,319 superTranscripts. The species-specific gene candidate\ntables contain the details of the overlaps with the reference\nannotations. Some superTranscripts overlap more than one gene, so\nI've retained the alternatives rather than assign a single gene name\nin every case. These are candidate gene associations, not confirmed\northologues. The original Necklace IDs have been retained throughout;\nthe lookup table links these to the IDs used for Liftoff.\n\nThe mapping and projection summaries are also in tables/, and the\noriginal superTranscript sequence and relevant annotations are in\nsource/.\n\n## Remaining annotation work\n\nThe functional annotation and statistical results still need to be\nadded. The integration script in scripts/ accepts the full Swiss-Prot\nand D. melanogaster-specific Trinotate reports separately.\n\nOne thing to check when adding the limma results: there are 142,593\nfeatures in the original graph annotation/count matrix and 121,757\nrows in TESTfeaturecounts.txt, the latter matching the number in\nthe manuscript. We still need to confirm the correspondence between\nblocks.gtf and the final limma feature IDs (including the row order\nused with make.unique()). The projected blocks.gtf track currently\nuses source row IDs, not the IDs generated by limma.\n\nIt would also be worth checking the group coding in the final R\nworkspace, as the comments at the beginning of Limma_alt_splice.R\nand the labels in the plotting section are reversed.\n\nI've kept the mapping scripts and the notes on the approaches we\ntried in https://github.com/MDSharma/simulans_choosiness. The copies\nin scripts/ are there for reference; the complete mapping pipeline\nand its original inputs are not included in this browser folder.\n\nIf needed, the original package files can be checked on Linux using\nsha256sum -c SHA256SUMS.\n\nMD\n",encoding="utf-8")
print("Wrote portable config.json and README.md",flush=True)
