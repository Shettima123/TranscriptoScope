# Release Notes

## 0.4.3 - 2026-06-11

- Made Windows package installation explicitly non-interactive by forcing binary package installs and suppressing source-package prompts.
- Added package-install progress messages for cleaner Windows runtime validation.

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
- Added built-in GO and KEGG pathway selection for yeast, human, and fruit fly.
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
- Pre-cached KEGG pathway mappings for yeast, human, and fruit fly.
- Improved Windows installer metadata and shortcut registration.
- Added safer uninstall path checks.

## 0.1.0

- Initial Windows-friendly Shiny app for raw count, normalized expression, and fold-change table analysis.
