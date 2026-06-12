# TranscriptoScope

A Windows-friendly Shiny app for differential gene expression analysis with DESeq2, normalized expression analysis, ORA enrichment, and ranked pathway analysis.

Created by Dr. Abubakar Abdulkadir, Southern University A and M.

This release can be installed for the current Windows user with a double-click installer. It lets a user upload raw counts, normalized expression matrices, or fold-change tables, inspect QC/results, run ORA or ranked pathway analysis, and download result files.

## Public Release

- Download the Windows zip from the GitHub Releases page.
- Installation instructions are in [INSTALL_WINDOWS.md](INSTALL_WINDOWS.md).
- Release notes are in [RELEASE_NOTES.md](RELEASE_NOTES.md).
- Citation metadata is in [CITATION.cff](CITATION.cff).
- License terms are in [LICENSE](LICENSE).
- Privacy notes are in [PRIVACY.md](PRIVACY.md).
- Support instructions are in [SUPPORT.md](SUPPORT.md).

If you use TranscriptoScope in research, cite the specific software release and cite the underlying statistical packages used in your analysis, including DESeq2 for differential expression and fgsea for ranked pathway analysis.

## What It Does

- Upload raw integer count matrices as `.csv`, `.tsv`, or `.txt`
- Upload normalized expression matrices with decimal values
- Upload fold-change tables with or without adjusted p-values
- Upload sample metadata as `.csv`, `.tsv`, or `.txt`, or let the app infer sample groups from names such as `condition_Rep1`
- Validate matching sample names between counts and metadata
- Select condition, reference group, comparison group, and optional batch column
- Run DESeq2 with `~ condition` or `~ batch + condition`
- Analyze normalized expression matrices with group means, log2 fold changes, Welch t-test p-values, and BH-adjusted p-values
- Visualize fold-change tables without requiring raw samples
- Export differential expression results and normalized counts
- Generate PCA, volcano, and MA plots
- Run over-representation analysis (ORA) from the current DEG result using built-in Ensembl GO gene sets, cached KEGG pathways, or custom gene sets
- Run Rosby's Lab-style ORA with separate up/down enrichment, protein-coding background, minimum overlap 2, default enrichment FDR 0.01, and top 5 terms per direction
- Run preranked pathway analysis with the full current DGE result, bundled GO or KEGG collections, normalized enrichment scores, FDR values, and leading-edge genes
- Plot top positive and negative pathways, a running enrichment curve, and a leading-edge expression heatmap when sample-level values are available
- Annotate DEG tables with bundled Ensembl mappings for yeast, human, and fruit fly
- Includes small bundled example data for testing

## Input Modes

### Read Counts Data

Use this for raw RNA-seq read counts. This is the true DESeq2 workflow.

The count matrix must contain raw integer counts. Do not use TPM, FPKM, RPKM, CPM, or already-normalized values.

The app also accepts converted count files with annotation columns such as `User_ID`, `ensembl_ID`, and `symbol`; nonnumeric annotation columns after the first gene ID column are ignored.

Example count matrix:

```csv
gene_id,ctrl_1,ctrl_2,treat_1,treat_2
GeneA,100,110,240,260
GeneB,55,61,52,59
```

Example metadata:

```csv
sample_id,condition,batch
ctrl_1,control,A
ctrl_2,control,B
treat_1,treated,A
treat_2,treated,B
```

The sample names in the count-matrix columns must match the `sample_id` values in metadata.

### Normalized Expression Data

Use this for decimal expression matrices such as normalized counts, TPM-like values, or other already-normalized sample-level measurements. This mode does not run DESeq2. It calculates group means, log2 fold change, Welch t-test p-values when replicates are available, and BH-adjusted p-values.

Example:

```csv
gene_id,S3_Galacose_Rep1,S3_Galacose_Rep2,S3_Glucose_Rep1,S3_Glucose_Rep2
YAL001C,1930.2,1906.1,2056.4,2076.8
YAL002W,2502.4,2727.4,1398.9,1618.2
```

If no metadata file is supplied, the app infers groups from sample names by removing suffixes such as `_Rep1`, `_Rep2`, or `_R1`.

### Fold Changes And Adjusted P-Values

Use this when another tool already produced a result table. Accepted adjusted p-value column names include `padj`, `FDR`, `qvalue`, and `adjusted_p_value`.

Example:

```csv
gene_id,log2FoldChange,pvalue,padj
GeneA,1.4,0.001,0.01
GeneB,-0.8,0.03,0.08
```

### Fold Changes Only

Use this for ranked fold-change lists with no p-values. The app can plot and export the ranking, but it cannot make significance calls.

## Enrichment / ORA

After running the analysis, open the Enrichment tab. The app automatically uses the current DEG result table.

By default, ORA uses a selected organism and gene set collection instead of requiring a gene set upload. `Auto from selected annotation` uses the organism selected in the Annotation panel.

The `Standard ORA` mode lets you choose up-regulated, down-regulated, or both DEG lists. It uses the current result genes as the enrichment universe and reports all tested terms sorted by adjusted p-value.

The `Rosby's Lab-style ORA` mode uses the lab's focused enrichment workflow:

- uses protein-coding genes as the background when a built-in annotation is selected
- intersects every gene set with that protein-coding background before calculating pathway size and enrichment
- tests upregulated and downregulated DEG genes separately
- uses minimum overlap 2
- uses FDR correction and defaults to enrichment FDR < 0.01
- displays up to 5 terms per direction
- exports lab-style enrichment columns

The app uses the selected built-in GO/KEGG source or custom gene set file, so results depend on the chosen collection and its release.

Supported built-in collections are:

- Ensembl BioMart GO mappings for yeast, human, and fruit fly, with optional filtering to biological process, molecular function, or cellular component
- KEGG pathways for yeast, human, and fruit fly, cached in `gene_sets/cache/`

The custom upload option is only for advanced use with your own gene sets. It is not for the count matrix. Custom gene set inputs can be:

- GMT files: `term_id`, `description`, then one or more gene IDs per line
- CSV/TSV long tables with columns such as `term_id`, `term_name`, and `gene_id`
- CSV/TSV wide tables with one row per gene set and genes across columns

The bundled `sample_data/gene_sets.csv` is only for testing the custom-upload parser. For real custom analysis, use gene sets that match the organism and gene IDs in your result table.

## Pathway Analysis

After running the DGE analysis, open the Pathway Analysis tab. The app automatically uses the current result table and ranks genes for preranked gene set enrichment analysis.

The default ranking is log2 fold change. Signed `-log10(p-value)` is also available when p-values are present. `Maximum gene FDR included` defaults to `1`, which keeps the full ranked list; lowering it filters genes before pathway analysis and should be done deliberately.

Choose the organism and a bundled GO or KEGG collection, then set the pathway size and FDR limits. The results include:

- enrichment score and normalized enrichment score
- raw and adjusted pathway p-values
- positive or negative enrichment direction
- leading-edge gene count and gene list
- top-pathway summary plot
- running enrichment plot for the selected pathway
- row-standardized leading-edge expression heatmap for raw-count or normalized-expression workflows

Fold-change table inputs can run pathway analysis, but they do not contain sample-level values, so the heatmap is unavailable. The Pathway Analysis tab uses the Bioconductor `fgsea` package.

## Built-In Ensembl Annotation

The Annotation panel can add gene symbol, description, biotype, chromosome, start, and end columns to the result table. The bundled annotations are:

- Yeast: `scerevisiae_gene_ensembl`, Saccharomyces cerevisiae R64-1-1, taxonomy `559292`
- Human: `hsapiens_gene_ensembl`, Homo sapiens GRCh38.p14, taxonomy `9606`
- Fruit fly: `dmelanogaster_gene_ensembl`, Drosophila melanogaster BDGP6.54, taxonomy `7227`

The app can match result IDs automatically against Ensembl/SGD IDs or gene symbols. Full annotations are included offline in the `annotations/` folder.

## Built-In GO And KEGG Gene Sets

The Enrichment tab includes offline GO term-to-gene mappings generated from Ensembl BioMart:

- Yeast: `scerevisiae_gene_ensembl`, Saccharomyces cerevisiae R64-1-1
- Human: `hsapiens_gene_ensembl`, Homo sapiens GRCh38.p14
- Fruit fly: `dmelanogaster_gene_ensembl`, Drosophila melanogaster BDGP6.54

Full GO mappings are included offline in the `gene_sets/` folder.

KEGG pathway mappings for yeast, human, and fruit fly are included in `gene_sets/cache/`. The app can refresh a missing cache from KEGG REST if needed.

## Install On Windows

1. Install R for Windows if it is not already installed.
2. Extract the app zip folder.
3. Double-click `Install_TranscriptoScope.bat`.
4. Launch the app from the Desktop or Start Menu shortcut.

The installer copies the app to `%LOCALAPPDATA%\TranscriptoScope`, checks the required R packages, and creates shortcuts. Use `Uninstall_TranscriptoScope.bat` or the Start Menu uninstall shortcut to remove it.

See `INSTALL_WINDOWS.md` for the short install guide.

The installed shortcut opens the app in an Edge or Chrome app-style window when available. The app is still powered by a local Shiny server on `127.0.0.1`; the quiet launcher hides the server window and opens the app interface directly.

## Run Locally On Windows Without Installing

1. Install R for Windows if it is not already installed.
2. Double-click `Install_Packages.bat`.
3. Double-click `Launch_TranscriptoScope.vbs`.

You can also run from PowerShell:

```powershell
Rscript scripts\install_packages.R
Rscript scripts\launch_app.R
```

## Developer Smoke Test

After dependencies are installed:

```powershell
Rscript scripts\smoke_test.R
```

The smoke test runs DESeq2, ORA, ranked pathway analysis, and pathway plot generation on the bundled example dataset without launching the Shiny interface.

## Packaging Direction

For early testers, zip this folder and ask users to run `Install_TranscriptoScope.bat`.

For a true installer, see `WINDOWS_PACKAGING.md`. The most practical path is:

1. Finish and test the Shiny app.
2. Freeze R/Bioconductor package versions with `renv`.
3. Create an installer with Inno Setup or RInno.
4. Code-sign the installer before public distribution.

## Files

- `app.R` - Shiny user interface and server logic
- `R/deseq_helpers.R` - reusable data validation and DESeq2 workflow helpers
- `sample_data/` - bundled example files for all input modes
- `sample_data/gene_sets.csv` - small demo gene set file for ORA testing
- `annotations/` - bundled Ensembl annotation CSVs and manifest
- `gene_sets/` - bundled Ensembl GO term-to-gene mappings, KEGG source manifest entries, and included KEGG cache
- `scripts/install_packages.R` - dependency installer
- `scripts/launch_app.R` - local app launcher
- `scripts/launch_app_window.ps1` - quiet Windows app-window launcher
- `scripts/download_ensembl_annotations.R` - regeneration script for bundled Ensembl annotations
- `scripts/download_ensembl_gene_sets.R` - regeneration script for bundled Ensembl GO gene sets
- `scripts/download_kegg_gene_sets.R` - optional script to pre-cache KEGG pathway mappings from KEGG REST
- `scripts/smoke_test.R` - command-line workflow test
- `scripts/build_release_zip.ps1` - creates the distributable Windows ZIP
- `Install_TranscriptoScope.bat` - Windows current-user installer
- `Launch_TranscriptoScope.vbs` - no-console Windows launcher used by shortcuts
- `Uninstall_TranscriptoScope.bat` - Windows current-user uninstaller
- `INSTALL_WINDOWS.md` - short user install guide
- `RELEASE_NOTES.md` - release summary
- `VERSION` - package version
- `WINDOWS_PACKAGING.md` - packaging and release guide
- `THIRD_PARTY_NOTES.md` - dependency/license notes
