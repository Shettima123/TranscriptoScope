# Resource Licensing And Redistribution

This document summarizes what is distributed with the public TranscriptoScope release and what is generated only on a user's local machine.

This is project documentation, not legal advice. Before institutional, commercial, or journal submission use, confirm the terms that apply to your specific distribution and audience.

| Resource type | Distributed in public ZIP? | Location | License/terms note |
| --- | --- | --- | --- |
| TranscriptoScope source code and project-authored documentation | Yes | Root folder, `R/`, `scripts/`, docs | MIT license |
| R and Bioconductor packages | No, installed by user or installer scripts | User R library | Each package keeps its own license and citation requirements |
| Ensembl-derived gene annotation CSVs | Yes | `annotations/` | Subject to Ensembl/EMBL-EBI terms and attribution guidance; source and timestamp are recorded in `annotations/manifest.csv` |
| Ensembl BioMart-derived GO mapping CSVs | Yes | `gene_sets/*_go.csv` | GO data products are CC BY 4.0; preserve GO attribution and release/source information |
| KEGG pathway mapping CSVs | No | `gene_sets/cache/` contains only `README.md` in the release | KEGG-derived cache files are generated only after user-initiated KEGG REST access and are not covered by MIT |
| User analysis outputs | No | Chosen by user during export | Owned/controlled by the user and subject to their institution's policies |

## KEGG Policy Used For This Release

The public release does not redistribute KEGG-derived pathway mapping CSV files. KEGG support is implemented as optional user-initiated access to KEGG REST. If the user selects a KEGG collection, the app downloads the selected organism's pathway mappings and stores them as a local cache on that user's machine.

Generated KEGG cache files should not be committed to the public repository, placed in release ZIPs, or redistributed unless the redistributor has confirmed permission or licensing for that use case.

The release build script removes `gene_sets/cache/*_kegg.csv` from the distributable ZIP as a safeguard.
