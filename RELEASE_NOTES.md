# Release Notes

## 0.4.5 - 2026-06-24

- Added `TF.Target.GTRD` as a built-in pathway database option for both enrichment and ranked pathway analysis.
- Added a Network Analysis tab with WGCNA module summaries, module-trait correlations, gene-module assignments, eigengenes, CSV downloads, and result-bundle exports.
- Added pathway cnetplot rendering, PNG export, cnetplot edge CSV export, and result-bundle outputs.
- Updated pathway cnetplot gene labels to prefer readable gene symbols over systematic gene identifiers when annotation is available.
- Added GO DAG plotting for over-represented GO terms, including ontology-derived ancestor relationships, significance-level coloring, PNG export, node and edge CSV exports, and bundle outputs.
- Added offline GO ontology relationship tables and a reproducible download script for rebuilding the DAG resources.

## 0.4.4 - 2026-06-21

- Added a top-level Refresh button that reloads the current Shiny session.
- Added a Significant Pathways CSV download for ranked pathway-analysis results below the selected pathway FDR cutoff.
- Added a Top Plotted Pathways CSV download that exports the exact pathway subset used by the Top Pathways plot.
- Added `pathway_significant_pathways.csv` and `pathway_top_plotted_pathways.csv` to the result bundle when pathway analysis has been run.
- Shared the top-pathway selection logic between plotting and CSV export to keep displayed and downloaded pathway subsets aligned.

## 0.4.3 - 2026-06-11

- Made Windows package installation explicitly non-interactive by forcing binary package installs and suppressing source-package prompts.
- Added package-install progress messages for cleaner Windows runtime validation.
- Updated public-release licensing notes to clarify that MIT covers TranscriptoScope source/documentation only and that KEGG-derived cache files are not bundled in the public ZIP.
- Added run-specific `analysis_report.md`, `reproduce_analysis.Rmd`, and `reproduce_analysis.R` files to the result bundle, with exact exported analysis inputs, matched gene-set tables for ORA/pathway reruns, and recorded UI settings.

## 0.4.2 - 2026-06-11

- Corrected Rosby's Lab-style ORA so each gene-set size is calculated after intersection with the declared protein-coding background.
- Prevented genes outside the analysis universe from changing enrichment p-values.
- Added a regression smoke test for out-of-background gene-set members.
- Fixed Windows package setup to create and use a per-user R library instead of attempting to install packages into Program Files.
- Improved uninstall cleanup by stopping running TranscriptoScope R processes from the install folder before removing files.

## 0.4.1 - 2026-06-10

- Improved Pathway Results table spacing and horizontal scrolling.
- Prevented compact headings from wrapping one character at a time.
- Replaced internal pathway result column names with readable display headings.

## 0.4.0 - 2026-06-10

- Added a Pathway Analysis tab that automatically uses the current DGE result.
- Added preranked GSEA with the Bioconductor fgsea package.
- Added built-in GO selection and optional KEGG pathway selection for yeast, human, and fruit fly.
- Added controls for ranking metric, pathway size, pathway FDR, optional gene FDR filtering, and absolute ranking.
- Added normalized enrichment scores, leading-edge genes, pathway summary plots, running enrichment plots, expression heatmaps, and CSV downloads.
- Added pathway outputs to the downloadable result bundle.

## 0.3.0 - 2026-06-10

- Renamed the app to TranscriptoScope.
- Added creator and affiliation branding: Dr. Abubakar Abdulkadir, Southern University A and M.
- Updated Windows install folder, shortcuts, uninstall entry, and release ZIP name for TranscriptoScope.
- Renamed the focused enrichment workflow to Rosby's Lab-style ORA.

## 0.2.1 - 2026-06-10

- Changed Desktop and Start Menu shortcuts to open a quiet app-style window instead of leaving a terminal window open.
- Added a top-level Start Menu shortcut in addition to the app folder shortcut for easier Windows search/start visibility.

## 0.2.0 - 2026-06-10

- Added the focused ORA workflow with protein-coding background, separate up/down enrichment, minimum overlap 2, default enrichment FDR 0.01, and top 5 displayed terms per direction.
- Added focused enrichment CSV exports.
- Included offline Ensembl annotations and GO mappings for yeast, human, and fruit fly.
- Added optional KEGG pathway support for yeast, human, and fruit fly through KEGG REST and a local user cache.
- Improved Windows installer metadata and shortcut registration.
- Added safer uninstall path checks.

## 0.1.0

- Initial Windows-friendly Shiny app for raw count, normalized expression, and fold-change table analysis.
