# Microsoft Store Listing Draft

## App Name

TranscriptoScope

## Short Description

Local Windows app for differential expression, enrichment, pathway, and WGCNA analysis of RNA-seq data.

## Full Description

TranscriptoScope is a local Windows application for RNA-seq differential gene-expression analysis, functional enrichment, ranked pathway analysis, and WGCNA module analysis. It provides a graphical workflow for preparing expression matrices, checking metadata alignment, running DESeq2-style count-based analysis, exploring normalized-expression inputs, visualizing differential-expression results, and exporting reproducible result bundles.

The app is designed for researchers, students, and laboratory teams who need a practical local workflow without writing every analysis step by hand. It supports raw read-count analysis with DESeq2, exploratory normalized-expression workflows, fold-change result inputs, built-in annotation resources for supported organisms, ORA enrichment, ranked pathway analysis, WGCNA module analysis, result plots, CSV exports, and reproducibility bundles containing result files and R code.

TranscriptoScope runs locally on the user's Windows computer. Uploaded datasets and analysis outputs remain on the local machine unless the user chooses to share them.

## Key Features

- Local RNA-seq differential-expression workflow.
- DESeq2-based raw count analysis.
- Input preflight checks for decimal matrices, metadata mismatch, replicate structure, and sample naming.
- Support for raw counts, normalized expression data, fold-change and adjusted p-value tables, and fold-change-only tables.
- Built-in annotation and GO gene-set resources for supported organisms.
- Optional KEGG retrieval through user-side cache when licensing and internet access allow.
- ORA enrichment and Rosby's Lab-style ORA.
- Ranked pathway analysis with fgsea.
- WGCNA module analysis for sample-level count or normalized-expression workflows.
- PCA, volcano, MA, regulation-summary, enrichment, and pathway plots.
- Downloadable result bundles with CSV outputs and reproducibility R/R Markdown code.

## Privacy Summary

TranscriptoScope runs locally. The app does not require users to upload data to a cloud service for analysis. Input files and result bundles remain on the user's computer. Optional online retrieval may occur only for external resources selected by the user, such as KEGG pathway resources.

## Category Suggestions

- Medical
- Education
- Productivity

If Microsoft requires a single best category, use Education for the first submission unless a more specific scientific/research category is available in Partner Center.

## Search Terms

RNA-seq, DESeq2, differential expression, bioinformatics, genomics, transcriptomics, pathway analysis, enrichment analysis, ORA, fgsea, WGCNA

## Support URL

https://github.com/Shettima123/TranscriptoScope/issues

## Website URL

https://github.com/Shettima123/TranscriptoScope

## Privacy Policy URL

https://github.com/Shettima123/TranscriptoScope/blob/main/PRIVACY.md

## Package Details

Package URL:

`https://tscope743176685.blob.core.windows.net/releases/TranscriptoScope_Windows_Store_v0.4.3.exe`

Architecture:

`x64`

Installer parameters:

`/VERYSILENT /SUPPRESSMSGBOXES /NORESTART`

App type:

`EXE`

## Notes For Certification

TranscriptoScope installs per user and does not require administrator privileges. The Store installer bundles the application, R runtime, and required R package library for offline installation. KEGG cache files are not bundled. Optional external resource downloads are user-triggered.
