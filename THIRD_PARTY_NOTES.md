# Third-Party Notes

This app is a wrapper around R and Bioconductor packages. Before public distribution, verify dependency licenses for the exact versions you ship.

Core dependencies:

- R
- Shiny
- ggplot2
- Bioconductor
- DESeq2
- fgsea

DESeq2 is distributed through Bioconductor and is licensed under LGPL. If you bundle R packages directly inside an installer, preserve license files and package citations.
fgsea is distributed through Bioconductor under the MIT license and provides the preranked pathway-analysis engine.

Suggested app citation text:

> Differential expression analysis is performed with DESeq2 from Bioconductor. Users should cite DESeq2 when reporting results.

> Preranked gene set enrichment analysis is performed with fgsea from Bioconductor. Users should cite fgsea when reporting pathway results.

## Ensembl Annotation

Bundled gene annotation CSVs and GO gene set CSVs are generated from Ensembl BioMart:

- `scerevisiae_gene_ensembl`
- `hsapiens_gene_ensembl`
- `dmelanogaster_gene_ensembl`

Before public distribution, verify the current Ensembl terms and attribution guidance for the exact annotation release you ship. Keep `annotations/manifest.csv` and `gene_sets/manifest.csv` with the package so users can see the source dataset, assembly, and download timestamp.

## KEGG Pathways

This release includes pre-cached KEGG pathway mappings for yeast, human, and fruit fly under `gene_sets/cache/`. If a cache is missing, the app can query the KEGG REST API for the selected organism and recreate it locally.

Before distributing a package that includes cached KEGG-derived mappings, verify KEGG access, redistribution, citation, and licensing requirements for your intended audience. KEGG's public Kyoto University/GenomeNet services request academic use, while commercial/web-service use may require licensing through Pathway Solutions.
