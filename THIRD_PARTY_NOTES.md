# Third-Party Notes

This app is a wrapper around R and Bioconductor packages. The TranscriptoScope MIT license applies to the app source code and project-authored documentation only. It does not relicense R, Bioconductor packages, annotation databases, GO data, KEGG data, or other third-party resources.

Core dependencies:

- R
- Shiny
- ggplot2
- Bioconductor
- DESeq2
- fgsea
- WGCNA
- impute
- preprocessCore

DESeq2 is distributed through Bioconductor and is licensed under LGPL. If you bundle R packages directly inside an installer, preserve license files and package citations.
fgsea is distributed through Bioconductor under the MIT license and provides the preranked pathway-analysis engine.
WGCNA is distributed through CRAN and provides the weighted gene co-expression network analysis engine. The WGCNA package imports Bioconductor packages including impute and preprocessCore.

Suggested app citation text:

> Differential expression analysis is performed with DESeq2 from Bioconductor. Users should cite DESeq2 when reporting results.

> Preranked gene set enrichment analysis is performed with fgsea from Bioconductor. Users should cite fgsea when reporting pathway results.

> Weighted gene co-expression network analysis is performed with WGCNA. Users should cite WGCNA when reporting module-network results.

## Ensembl Annotation

Bundled gene annotation CSVs and GO gene set CSVs are generated from Ensembl BioMart:

- `scerevisiae_gene_ensembl`
- `hsapiens_gene_ensembl`
- `dmelanogaster_gene_ensembl`

Keep `annotations/manifest.csv` and `gene_sets/manifest.csv` with the package so users can see the source dataset, assembly, and download timestamp. Ensembl data and BioMart-derived exports remain subject to Ensembl/EMBL-EBI terms and attribution guidance.

GO term-to-gene mappings are generated from BioMart GO annotation fields. Gene Ontology Consortium data products are licensed under CC BY 4.0, and public redistribution should preserve GO attribution, release/version information where available, and license notice.

GO DAG plots use compact term and relationship CSV files generated from the Gene Ontology `go-basic.obo` file. These files are distributed under the same GO data product terms and can be regenerated with `scripts/download_go_ontology.R`.

## MSigDB TF.Target.GTRD

The bundled human `TF.Target.GTRD` gene sets are generated from the MSigDB C3:TFT GTRD symbol GMT file. MSigDB and GTRD resources remain third-party biological data and are not covered by the TranscriptoScope MIT license. Preserve the source URL and release recorded in `gene_sets/manifest.csv`, and confirm that your intended redistribution or publication use follows MSigDB and GTRD terms.

## KEGG Pathways

This public release does not bundle KEGG pathway mapping CSV files. The `gene_sets/cache/` directory contains only a README in the release package.

When a user selects KEGG, the app can query KEGG REST for the selected organism and cache the resulting term-to-gene table locally on that user's machine. These locally generated cache files are KEGG-derived data and are not covered by the TranscriptoScope MIT license. Do not redistribute generated KEGG cache files unless you have confirmed that redistribution is permitted for your use case.

KEGG's public legal page states that academic users may freely use the KEGG website, while non-academic use requires a commercial license. Pathway Solutions provides KEGG commercial licensing, including licenses that can cover web/API access. Users and redistributors are responsible for following the KEGG terms that apply to their institution and use case.
