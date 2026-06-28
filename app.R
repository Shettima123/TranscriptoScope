if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The shiny package is required. Run scripts/install_packages.R, then launch again.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required. Run scripts/install_packages.R, then launch again.", call. = FALSE)
}

library(shiny)

source(file.path("R", "deseq_helpers.R"), local = TRUE)

app_version <- tryCatch(
  trimws(readLines("VERSION", warn = FALSE)[1]),
  error = function(err) "0.4.6"
)
if (!nzchar(app_version) || is.na(app_version)) {
  app_version <- "0.4.6"
}

matrix_modes <- c("counts", "normalized")
annotation_dir <- file.path("annotations")
annotation_manifest <- tryCatch(read_annotation_manifest(annotation_dir), error = function(err) data.frame())
annotation_choices <- c("None" = "none")
if (nrow(annotation_manifest) > 0) {
  annotation_choices <- c(annotation_choices, stats::setNames(annotation_manifest$key, annotation_manifest$label))
}
default_annotation <- if ("yeast_scerevisiae" %in% annotation_choices) "yeast_scerevisiae" else "none"
gene_set_dir <- file.path("gene_sets")
gene_set_manifest <- tryCatch(read_gene_set_manifest(gene_set_dir), error = function(err) data.frame())
gene_set_choices <- c("Auto from selected annotation" = "auto")
if (nrow(gene_set_manifest) > 0) {
  organism_rows <- gene_set_manifest[!duplicated(gene_set_manifest$annotation_key), , drop = FALSE]
  organism_labels <- organism_rows$scientific_name
  if (nrow(annotation_manifest) > 0) {
    matched_labels <- annotation_manifest$label[match(organism_rows$annotation_key, annotation_manifest$key)]
    organism_labels[!is.na(matched_labels) & matched_labels != ""] <- matched_labels[!is.na(matched_labels) & matched_labels != ""]
  }
  gene_set_choices <- c(gene_set_choices, stats::setNames(organism_rows$annotation_key, organism_labels))
}
gene_set_choices <- c(gene_set_choices, "Custom gene set file" = "custom")
pathway_gene_set_choices <- gene_set_choices[gene_set_choices != "custom"]
gene_set_collection_choices <- c(
  "GO: all terms" = "go_all",
  "GO: biological process" = "go_biological_process",
  "GO: molecular function" = "go_molecular_function",
  "GO: cellular component" = "go_cellular_component",
  "KEGG pathways" = "kegg_pathway",
  "TF.Target.GTRD" = "tf_target_gtrd"
)
enrichment_mode_choices <- c(
  "Standard ORA" = "standard",
  "Rosby's Lab-style ORA" = "rosbys_lab"
)

mode_label <- function(mode) {
  switch(
    mode,
    counts = "Read counts data (DESeq2)",
    normalized = "Normalized expression data",
    fc_padj = "Fold changes and adjusted p-values",
    fc_only = "Fold changes only",
    "Input data"
  )
}

example_file_for_mode <- function(mode) {
  switch(
    mode,
    counts = "counts.csv",
    normalized = "normalized_expression.csv",
    fc_padj = "fold_changes_padj.csv",
    fc_only = "fold_changes_only.csv",
    "counts.csv"
  )
}

example_download_specs <- list(
  list(id = "download_example_counts", label = "Example Counts CSV", file = "counts.csv"),
  list(id = "download_example_metadata", label = "Example Metadata CSV", file = "metadata.csv"),
  list(id = "download_example_normalized", label = "Example Normalized CSV", file = "normalized_expression.csv"),
  list(id = "download_example_fc_padj", label = "Example FC + padj CSV", file = "fold_changes_padj.csv"),
  list(id = "download_example_fc_only", label = "Example FC Only CSV", file = "fold_changes_only.csv"),
  list(id = "download_example_gene_sets", label = "Example Gene Sets CSV", file = "gene_sets.csv"),
  list(id = "download_example_bundle", label = "All Example Inputs ZIP", file = NA_character_)
)

gene_set_collection_info <- function(collection_choice, default = "go_all") {
  collection_choice <- if (is.null(collection_choice) || !nzchar(collection_choice)) default else collection_choice
  switch(
    collection_choice,
    go_biological_process = list(collection = "go", domain = "biological_process"),
    go_molecular_function = list(collection = "go", domain = "molecular_function"),
    go_cellular_component = list(collection = "go", domain = "cellular_component"),
    kegg_pathway = list(collection = "kegg", domain = "kegg_pathway"),
    tf_target_gtrd = list(collection = "gtrd", domain = "tf_target_gtrd"),
    list(collection = "go", domain = "all")
  )
}

diagnostic_details <- function(message, input_name = NULL, input_type = NULL, context = NULL, probable_cause = NULL) {
  message <- trimws(as.character(message)[1])
  probable_cause <- trimws(as.character(probable_cause)[1])
  if (is.na(probable_cause)) {
    probable_cause <- ""
  }
  if (is.na(message) || !nzchar(message)) {
    message <- if (nzchar(probable_cause)) probable_cause else "TranscriptoScope could not identify the exact error message."
  }
  input_name <- if (is.null(input_name) || is.na(input_name) || !nzchar(input_name)) "" else input_name
  input_label <- mode_label(if (is.null(input_type)) "" else input_type)
  text <- tolower(paste(message, probable_cause, input_name, input_label, context, collapse = " "))

  default_fixes <- c(
    "Confirm that the selected input data type matches the file you uploaded.",
    "Use the Preflight tab to check gene count, sample count, metadata rows, and warnings before running analysis.",
    "If the file came from Excel, export the active worksheet as CSV or tab-delimited text before uploading."
  )

  make_details <- function(title, fixes) {
    list(
      title = title,
      message = message,
      probable_cause = if (nzchar(probable_cause) && !grepl(probable_cause, message, fixed = TRUE)) probable_cause else "",
      fixes = unique(c(fixes, default_fixes))
    )
  }

  if (grepl("excel|\\.xlsx|\\.xls|embedded null|eof within quoted string", text)) {
    return(make_details(
      "Excel workbook detected",
      c(
        "Open the workbook in Excel and save the worksheet you want as CSV or tab-delimited text.",
        "Upload the exported CSV/TSV file as the data file.",
        "If the worksheet contains decimal CPM/TPM/FPKM/logCPM values, choose Normalized expression data instead of Read counts data (DESeq2)."
      )
    ))
  }

  if (grepl("must contain at least two columns|needs a gene id column plus one or more", text)) {
    return(make_details(
      "The file does not look like an analysis table",
      c(
        "Use a header row with gene IDs in the first column.",
        "Add at least one numeric sample column for matrix inputs, or a fold-change column for fold-change modes.",
        "For CSV files, use commas; for TSV/TXT files, use tabs."
      )
    ))
  }

  if (grepl("decimal values", text)) {
    return(make_details(
      "Raw-count mode received decimal values",
      c(
        "Choose Normalized expression data if the values are CPM, TPM, FPKM, RPKM, normalized counts, or log expression.",
        "Use Read counts data (DESeq2) only for raw non-negative integer counts from tools such as featureCounts or HTSeq.",
        "Do not round normalized expression values to force DESeq2; use the original raw count matrix when available."
      )
    ))
  }

  if (grepl("missing or non-numeric|blank or non-numeric|empty values|contains missing|non-finite|cannot be used", text)) {
    return(make_details(
      "The file contains blank or non-numeric values",
      c(
        "Find and fill or remove blank cells in numeric sample columns.",
        "Remove text such as comments, percent signs, or formulas that produce blanks from numeric columns.",
        "If only one gene/sample has the problem, remove that row or correct the cell and upload again."
      )
    ))
  }

  if (grepl("gene ids must be unique|first duplicate|sample names must be unique|metadata sample ids must be unique", text)) {
    return(make_details(
      "Duplicate identifiers were found",
      c(
        "Make gene IDs unique before upload, or collapse duplicate gene rows outside the app.",
        "Make sample column names unique.",
        "For metadata, each sample_id must appear only once."
      )
    ))
  }

  if (grepl("metadata|sample id|sample\\(s\\)", text)) {
    return(make_details(
      "Metadata does not match the data matrix",
      c(
        "Make the metadata sample_id values exactly match the sample column names in the data file.",
        "Check spelling, spaces, capitalization, punctuation, and hidden trailing spaces.",
        "If you do not need custom metadata, remove the metadata file and let the app infer groups from sample names."
      )
    ))
  }

  if (grepl("condition column|reference group|comparison group|at least two groups|replicates", text)) {
    return(make_details(
      "The selected comparison cannot be formed",
      c(
        "Choose a condition column with at least two groups.",
        "Choose different reference and comparison groups.",
        "For statistical testing, include biological replicates in each group whenever possible."
      )
    ))
  }

  if (grepl("confounded|cannot be fit|simpler design|interaction factor|design factor", text)) {
    return(make_details(
      "The selected DESeq2 design is not estimable",
      c(
        "Remove adjustment factors that perfectly match the condition groups.",
        "Use a simpler design when a batch, donor, or interaction factor is confounded with condition.",
        "Check the metadata table to confirm every factor level has samples in the needed condition groups."
      )
    ))
  }

  if (grepl("no genes pass|fewer than two genes pass|minimum total count|minimum cpm|zero total counts", text)) {
    return(make_details(
      "Filtering removed too much data",
      c(
        "Lower the minimum total count or CPM threshold.",
        "Check that sample columns contain expression/count values and not metadata.",
        "Remove samples with zero total counts before uploading."
      )
    ))
  }

  if (grepl("fold-change|adjusted p-value|padj|fdr|qvalue", text)) {
    return(make_details(
      "Fold-change table columns could not be recognized",
      c(
        "Include a gene ID column and a numeric fold-change column.",
        "For Fold changes and adjusted p-values mode, include a padj, FDR, qvalue, or adjusted_p_value column.",
        "Use log2 fold changes unless you intentionally select linear fold-change scale."
      )
    ))
  }

  if (grepl("gene sets overlapped|pathways overlap|ranked genes|size filtering|selected deg genes|enrichment", text)) {
    return(make_details(
      "Enrichment/pathway inputs do not overlap enough genes",
      c(
        "Confirm the selected organism matches the gene IDs in the result table.",
        "Try the other gene-ID matching mode or turn off case-sensitive matching.",
        "Relax DEG, pathway-size, or gene-set-size filters and rerun."
      )
    ))
  }

  if (grepl("kegg|download|internet|rest", text)) {
    return(make_details(
      "The selected online resource could not be loaded",
      c(
        "Confirm the computer has internet access if KEGG is selected.",
        "Try GO or TF.Target.GTRD if you need a fully bundled/offline source.",
        "Restart the app and retry after the network connection is stable."
      )
    ))
  }

  if (grepl("package is not installed|required", text)) {
    return(make_details(
      "A required R package is missing",
      c(
        "Run Install_Packages.bat from the TranscriptoScope folder.",
        "Restart TranscriptoScope after package installation finishes.",
        "If installation fails, check your internet connection and R package library permissions."
      )
    ))
  }

  if (grepl("wgcna", text)) {
    return(make_details(
      "WGCNA requirements were not met",
      c(
        "Use count or normalized-expression input, not fold-change-only input.",
        "Use at least four samples and enough variable genes after filtering.",
        "Lower the WGCNA variable-gene filter or use a larger sample-level expression matrix."
      )
    ))
  }

  make_details("TranscriptoScope could not continue", character())
}

diagnose_uploaded_table <- function(input_path = NULL, input_name = NULL, input_type = NULL) {
  if (is.null(input_path) || is.na(input_path) || !nzchar(input_path) || !file.exists(input_path)) {
    return(NULL)
  }
  input_name <- if (is.null(input_name) || is.na(input_name) || !nzchar(input_name)) basename(input_path) else input_name
  lower_name <- tolower(input_name)
  if (grepl("\\.(xlsx|xls)$", lower_name)) {
    return("The uploaded file is an Excel workbook. Save the worksheet as CSV or TSV, then upload that exported text file.")
  }

  data <- tryCatch(
    read_dge_table(input_path, input_name),
    error = function(err) {
      sprintf("The uploaded file could not be read as a table: %s", conditionMessage(err))
    }
  )
  if (is.character(data)) {
    return(data)
  }
  if (ncol(data) < 2) {
    return("The uploaded table has fewer than two columns. Use gene IDs in the first column and sample/result values in additional columns.")
  }

  gene_ids <- trimws(as.character(data[[1]]))
  if (any(is.na(gene_ids) | gene_ids == "")) {
    row_index <- which(is.na(gene_ids) | gene_ids == "")[1] + 1L
    return(sprintf("The first gene-ID column has a blank value at data row %s.", row_index))
  }
  if (anyDuplicated(gene_ids)) {
    duplicated_gene <- gene_ids[duplicated(gene_ids)][1]
    return(sprintf("The uploaded table has duplicate gene IDs. First duplicate: %s.", duplicated_gene))
  }

  value_columns <- data[-1]
  numeric_values <- lapply(value_columns, function(column) {
    values <- trimws(as.character(column))
    values[is.na(values)] <- ""
    suppressWarnings(as.numeric(gsub(",", "", values, fixed = TRUE)))
  })
  numeric_column_has_values <- vapply(numeric_values, function(column) any(!is.na(column)), logical(1))
  if (!any(numeric_column_has_values)) {
    return("No numeric value columns were detected after the gene-ID column.")
  }

  for (column_index in which(numeric_column_has_values)) {
    raw_values <- trimws(as.character(value_columns[[column_index]]))
    raw_values[is.na(raw_values)] <- ""
    numeric_column <- numeric_values[[column_index]]
    bad_rows <- which(raw_values == "" | is.na(numeric_column))
    if (length(bad_rows) > 0) {
      row_index <- bad_rows[1]
      return(sprintf(
        "The uploaded table contains a blank or non-numeric value at row %s, gene %s, column %s.",
        row_index + 1L,
        gene_ids[row_index],
        names(value_columns)[column_index]
      ))
    }
  }

  if (identical(input_type, "counts")) {
    for (column_index in which(numeric_column_has_values)) {
      numeric_column <- numeric_values[[column_index]]
      decimal_rows <- which(!is.na(numeric_column) & abs(numeric_column - round(numeric_column)) > 1e-6)
      if (length(decimal_rows) > 0) {
        row_index <- decimal_rows[1]
        return(sprintf(
          "Read counts mode found decimal values, starting at row %s, gene %s, column %s. These look like normalized expression values rather than raw integer counts.",
          row_index + 1L,
          gene_ids[row_index],
          names(value_columns)[column_index]
        ))
      }
    }
  }

  NULL
}

diagnostic_alert <- function(message, severity = "warning", input_name = NULL, input_path = NULL, input_type = NULL, context = NULL) {
  probable_cause <- diagnose_uploaded_table(
    input_path = input_path,
    input_name = input_name,
    input_type = input_type
  )
  details <- diagnostic_details(
    message = message,
    input_name = input_name,
    input_type = input_type,
    context = context,
    probable_cause = probable_cause
  )
  div(
    class = paste("alert-block", severity, "diagnostic-block"),
    tags$strong(details$title),
    tags$p(class = "diagnostic-message", details$message),
    if (nzchar(details$probable_cause)) {
      tags$p(class = "diagnostic-message", tags$strong("Most likely cause: "), details$probable_cause)
    },
    tags$div(
      class = "diagnostic-fixes",
      tags$span("How to address it"),
      tags$ul(lapply(details$fixes, tags$li))
    )
  )
}

app_css <- "
:root {
  --ink: #162033;
  --muted: #526071;
  --line: #D8DEE8;
  --panel: #FFFFFF;
  --surface: #F5F7FA;
  --accent: #0F766E;
  --accent-dark: #115E59;
  --danger: #B42318;
  --warning: #8A5A00;
}

body {
  background: var(--surface);
  color: var(--ink);
  font-family: 'Segoe UI', system-ui, -apple-system, BlinkMacSystemFont, sans-serif;
}

.app-shell {
  max-width: 1480px;
  margin: 0 auto;
  padding: 22px;
}

.app-header {
  display: flex;
  align-items: flex-end;
  justify-content: space-between;
  gap: 18px;
  padding-bottom: 16px;
  border-bottom: 1px solid var(--line);
  margin-bottom: 18px;
}

.app-title h1 {
  font-size: 28px;
  line-height: 1.15;
  margin: 0;
  font-weight: 700;
  letter-spacing: 0;
}

.title-row {
  display: flex;
  align-items: baseline;
  gap: 10px;
  flex-wrap: wrap;
}

.version-badge {
  display: inline-flex;
  align-items: center;
  border: 1px solid #99D0C9;
  border-radius: 999px;
  background: #E6F5F2;
  color: var(--accent-dark);
  font-size: 12px;
  font-weight: 700;
  line-height: 1;
  padding: 5px 9px;
}

.app-title p {
  margin: 6px 0 0;
  color: var(--muted);
  font-size: 14px;
}

.creator-line {
  color: var(--ink);
  font-weight: 600;
}

.status-line {
  font-size: 13px;
  color: var(--muted);
  white-space: nowrap;
}

.header-actions {
  display: flex;
  flex-direction: column;
  align-items: flex-end;
  gap: 8px;
}

.refresh-button {
  border-color: #99D0C9;
  color: var(--accent-dark);
  background: #F4FBF9;
  font-weight: 700;
  min-width: 92px;
}

.refresh-button:hover,
.refresh-button:focus {
  border-color: var(--accent);
  color: var(--accent-dark);
  background: #E6F5F2;
}

.layout-grid {
  display: grid;
  grid-template-columns: minmax(280px, 360px) minmax(0, 1fr);
  gap: 18px;
  align-items: start;
}

.layout-grid > *,
.sidebar,
.main {
  min-width: 0;
}

.panel {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 16px;
  box-shadow: 0 10px 24px rgba(22, 32, 51, 0.06);
  max-width: 100%;
  overflow-x: auto;
}

.panel + .panel {
  margin-top: 14px;
}

.panel h2,
.panel h3 {
  margin-top: 0;
  letter-spacing: 0;
}

.panel h2 {
  font-size: 18px;
  margin-bottom: 14px;
}

.panel h3 {
  font-size: 15px;
  margin-bottom: 10px;
}

.form-group {
  margin-bottom: 13px;
}

.btn-primary {
  background: var(--accent);
  border-color: var(--accent-dark);
}

.btn-primary:hover,
.btn-primary:focus {
  background: var(--accent-dark);
  border-color: var(--accent-dark);
}

.run-button {
  width: 100%;
  font-weight: 700;
}

.soft-note {
  color: var(--muted);
  font-size: 13px;
  line-height: 1.45;
}

.alert-block {
  border-left: 4px solid var(--accent);
  background: #ECFDF5;
  padding: 12px 14px;
  margin-bottom: 14px;
}

.alert-block.warning {
  border-left-color: var(--warning);
  background: #FFF7ED;
}

.alert-block.error {
  border-left-color: var(--danger);
  background: #FEF2F2;
}

.diagnostic-block strong {
  display: block;
  margin-bottom: 6px;
}

.diagnostic-message {
  margin: 0 0 8px;
}

.diagnostic-fixes span {
  display: block;
  font-weight: 700;
  margin-bottom: 4px;
}

.diagnostic-fixes ul {
  margin: 0;
  padding-left: 20px;
}

.diagnostic-fixes li {
  margin-bottom: 3px;
}

.metric-grid {
  display: grid;
  grid-template-columns: repeat(4, minmax(120px, 1fr));
  gap: 10px;
  margin-bottom: 14px;
}

.metric {
  border: 1px solid var(--line);
  border-radius: 8px;
  padding: 12px;
  min-height: 78px;
  background: #FBFCFE;
}

.metric span {
  display: block;
  color: var(--muted);
  font-size: 12px;
}

.metric strong {
  display: block;
  font-size: 22px;
  line-height: 1.2;
  margin-top: 6px;
}

.tab-content {
  padding-top: 14px;
}

table {
  font-size: 13px;
  max-width: 100%;
}

.nav-tabs {
  display: flex;
  flex-wrap: wrap;
}

.alert-block,
.soft-note,
td {
  overflow-wrap: anywhere;
}

th {
  overflow-wrap: normal;
  word-break: normal;
}

#pathway_results_table {
  width: 100%;
  overflow-x: auto;
}

#pathway_results_table table {
  width: 1180px !important;
  min-width: 1180px;
  max-width: none;
  table-layout: fixed;
}

#pathway_results_table th {
  white-space: normal;
  vertical-align: bottom;
  line-height: 1.25;
}

#pathway_results_table td {
  vertical-align: top;
  line-height: 1.4;
}

#pathway_results_table th:first-child,
#pathway_results_table td:first-child {
  width: 220px;
}

#pathway_results_table th:last-child,
#pathway_results_table td:last-child {
  width: 330px;
}

.download-grid {
  display: grid;
  grid-template-columns: repeat(3, minmax(180px, 1fr));
  gap: 12px;
}

.download-grid .btn {
  width: 100%;
}

@media (max-width: 920px) {
  .app-shell {
    padding: 14px;
  }

  .app-header {
    display: block;
  }

  .status-line {
    margin-top: 8px;
    white-space: normal;
  }

  .header-actions {
    align-items: flex-start;
    margin-top: 12px;
  }

  .layout-grid {
    grid-template-columns: 1fr;
  }

  .metric-grid,
  .download-grid {
    grid-template-columns: 1fr;
  }
}
"

ui <- fluidPage(
  tags$head(
    tags$title(sprintf("TranscriptoScope v%s", app_version)),
    tags$link(rel = "icon", type = "image/x-icon", href = "favicon.ico"),
    tags$style(HTML(app_css)),
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('transcriptoscope-refresh', function(message) {
        window.location.replace(window.location.pathname + window.location.search);
      });
      document.addEventListener('click', function(event) {
        var button = event.target.closest ? event.target.closest('#refresh_app') : null;
        if (!button) {
          return;
        }
        window.setTimeout(function() {
          window.location.replace(window.location.pathname + window.location.search);
        }, 800);
      });"
    ))
  ),
  div(
    class = "app-shell",
    div(
      class = "app-header",
      div(
        class = "app-title",
        div(
          class = "title-row",
          h1("TranscriptoScope"),
          span(class = "version-badge", sprintf("v%s", app_version))
        ),
        p("Windows-friendly transcriptomics analysis for raw counts, expression matrices, fold-change tables, enrichment, pathways, and WGCNA"),
        p(class = "creator-line", "Created by Dr. Abubakar Abdulkadir | Dr. Rosby's Lab, Southern University A and M")
      ),
      div(
        class = "header-actions",
        actionButton("refresh_app", "Reset / Refresh", class = "refresh-button", title = "Clear errors, reset the app state, and reload"),
        div(class = "status-line", textOutput("runtime_status", inline = TRUE))
      )
    ),
    div(
      class = "layout-grid",
      div(
        class = "sidebar",
        div(
          class = "panel",
          h2("Inputs"),
          selectInput(
            "input_type",
            "Input data type",
            choices = c(
              "Read counts data (DESeq2)" = "counts",
              "Normalized expression data" = "normalized",
              "Fold changes & adjusted p-values" = "fc_padj",
              "Fold changes only" = "fc_only"
            ),
            selected = "counts",
            selectize = FALSE
          ),
          checkboxInput("use_example", "Use bundled example data", value = FALSE),
          fileInput(
            "count_file",
            "Data file",
            accept = c(".csv", ".tsv", ".txt")
          ),
          div(
            class = "soft-note",
            "Use CSV, TSV, or TXT tables. If your data are in Excel, save the worksheet as CSV/TSV first. Decimal CPM/TPM/FPKM/logCPM values belong in Normalized expression data mode."
          ),
          uiOutput("metadata_upload"),
          uiOutput("input_status")
        ),
        div(
          class = "panel",
          h2("Annotation"),
          selectInput(
            "annotation_species",
            "Built-in Ensembl annotation",
            choices = annotation_choices,
            selected = default_annotation,
            selectize = FALSE
          ),
          selectInput(
            "annotation_match",
            "Match result gene IDs by",
            choices = c(
              "Auto-detect" = "auto",
              "Ensembl/SGD gene ID" = "ensembl",
              "Gene symbol" = "symbol"
            ),
            selected = "auto",
            selectize = FALSE
          ),
          uiOutput("annotation_status")
        ),
        div(
          class = "panel",
          h2("Analysis"),
          uiOutput("analysis_controls")
        )
      ),
      div(
        class = "main",
        uiOutput("analysis_status"),
        tabsetPanel(
          id = "main_tabs",
          tabPanel(
            "Preprocess",
            div(
              class = "panel",
              h2("Preprocess Counts"),
              uiOutput("preprocess_status"),
              div(
                class = "download-grid",
                numericInput("preprocess_min_cpm", "Minimum CPM", value = 0.5, min = 0, step = 0.1),
                numericInput("preprocess_min_samples", "Minimum samples", value = 1, min = 1, step = 1),
                numericInput("preprocess_pseudocount", "Transform pseudocount", value = 4, min = 0, step = 1)
              ),
              uiOutput("preprocess_metrics"),
              tableOutput("preprocess_sample_table"),
              div(
                class = "download-grid",
                downloadButton("download_preprocessed_counts", "Filtered Counts"),
                downloadButton("download_preprocessed_expression", "log2(CPM + c)"),
                downloadButton("download_preprocessed_metadata", "Metadata CSV")
              )
            )
          ),
          tabPanel(
            "Preflight",
            div(
              class = "panel",
              h2("Preflight"),
              uiOutput("preflight_metrics"),
              tableOutput("preflight_table"),
              uiOutput("preflight_warnings")
            )
          ),
          tabPanel(
            "Results",
            div(
              class = "panel",
              h2("Results"),
              tableOutput("results_table")
            )
          ),
          tabPanel(
            "Enrichment",
            div(
              class = "panel",
              h2("ORA Enrichment"),
              uiOutput("enrichment_status"),
              div(
                class = "download-grid",
                selectInput(
                  "enrichment_mode",
                  "Enrichment method",
                  choices = enrichment_mode_choices,
                  selected = "standard",
                  selectize = FALSE
                ),
                selectInput(
                  "gene_set_source",
                  "Gene set organism",
                  choices = gene_set_choices,
                  selected = "auto",
                  selectize = FALSE
                ),
                selectInput(
                  "go_domain",
                  "Gene set collection",
                  choices = gene_set_collection_choices,
                  selected = "go_all",
                  selectize = FALSE
                ),
                selectInput(
                  "ora_gene_list",
                  "DEG list",
                  choices = c(
                    "Up regulated" = "up",
                    "Down regulated" = "down",
                    "Up + Down regulated" = "both"
                  ),
                  selected = "both",
                  selectize = FALSE
                )
              ),
              div(
                class = "download-grid",
                numericInput("ora_min_set_size", "Minimum gene set size", value = 3, min = 1, step = 1),
                numericInput("ora_max_set_size", "Maximum gene set size", value = 2000, min = 2, step = 50),
                numericInput("ora_padj_cutoff", "Enrichment FDR cutoff", value = 0.1, min = 0.001, max = 1, step = 0.01)
              ),
              div(
                class = "download-grid",
                numericInput("ora_top_n", "Top terms to plot", value = 20, min = 5, max = 100, step = 5),
                checkboxInput("ora_case_sensitive", "Case-sensitive custom gene matching", value = FALSE)
              ),
              uiOutput("enrichment_mode_note"),
              uiOutput("custom_gene_set_upload"),
              tableOutput("enrichment_summary")
            ),
            div(
              class = "panel",
              h2("Top Enriched Gene Sets"),
              plotOutput("enrichment_plot", height = "430px")
            ),
            div(
              class = "panel",
              h2("GO DAG Plot"),
              uiOutput("go_dag_status"),
              div(
                class = "download-grid",
                numericInput("go_dag_top_n", "Top GO terms in DAG", value = 12, min = 1, max = 60, step = 1),
                numericInput("go_dag_ancestor_depth", "Ancestor levels", value = 3, min = 0, max = 8, step = 1),
                numericInput("go_dag_padj_cutoff", "DAG FDR cutoff", value = 0.1, min = 0.001, max = 1, step = 0.01),
                downloadButton("download_go_dag_plot", "GO DAG PNG"),
                downloadButton("download_go_dag_nodes", "GO DAG Nodes CSV"),
                downloadButton("download_go_dag_edges", "GO DAG Edges CSV")
              ),
              plotOutput("go_dag_plot", height = "720px")
            ),
            div(
              class = "panel",
              h2("Enrichment Results"),
              tableOutput("enrichment_table"),
              div(
                class = "download-grid",
                downloadButton("download_enrichment", "Enrichment CSV"),
                downloadButton("download_selected_ora_genes", "Selected Genes CSV")
              )
            )
          ),
          tabPanel(
            "Pathway Analysis",
            div(
              class = "panel",
              h2("Ranked Pathway Analysis"),
              uiOutput("pathway_status"),
              div(
                class = "download-grid",
                selectInput(
                  "pathway_gene_set_source",
                  "Gene set organism",
                  choices = pathway_gene_set_choices,
                  selected = "auto",
                  selectize = FALSE
                ),
                selectInput(
                  "pathway_collection",
                  "Pathway collection",
                  choices = gene_set_collection_choices,
                  selected = "kegg_pathway",
                  selectize = FALSE
                ),
                selectInput(
                  "pathway_rank_metric",
                  "Gene ranking metric",
                  choices = c(
                    "Log2 fold change" = "log2fc",
                    "Signed -log10 p-value" = "signed_p"
                  ),
                  selected = "log2fc",
                  selectize = FALSE
                )
              ),
              div(
                class = "download-grid",
                numericInput("pathway_min_set_size", "Minimum pathway size", value = 3, min = 1, step = 1),
                numericInput("pathway_max_set_size", "Maximum pathway size", value = 500, min = 2, step = 25),
                numericInput("pathway_padj_cutoff", "Pathway FDR cutoff", value = 0.05, min = 0.001, max = 1, step = 0.01)
              ),
              div(
                class = "download-grid",
                numericInput("pathway_gene_padj_cutoff", "Maximum gene FDR included", value = 1, min = 0.001, max = 1, step = 0.05),
                numericInput("pathway_top_n", "Top pathways to plot", value = 20, min = 5, max = 100, step = 5),
                div(
                  checkboxInput("pathway_absolute_ranking", "Use absolute ranking scores", value = FALSE),
                  checkboxInput("pathway_show_ids", "Show pathway IDs", value = FALSE)
                )
              ),
              actionButton(
                "run_pathway_analysis",
                "Run Pathway Analysis",
                class = "btn-primary run-button"
              ),
              tableOutput("pathway_summary")
            ),
            div(
              class = "panel",
              h2("Top Pathways"),
              plotOutput("pathway_summary_plot", height = "520px")
            ),
            div(
              class = "panel",
              h2("Pathway Cnetplot"),
              div(
                class = "download-grid",
                numericInput("pathway_cnet_top_n", "Top pathways in cnetplot", value = 5, min = 1, max = 20, step = 1),
                numericInput("pathway_cnet_max_genes", "Max genes per pathway", value = 15, min = 3, max = 100, step = 1),
                downloadButton("download_pathway_cnetplot", "Pathway Cnetplot PNG"),
                downloadButton("download_pathway_cnetplot_tiff", "Pathway Cnetplot TIFF"),
                downloadButton("download_pathway_cnet_edges", "Cnetplot Edges CSV")
              ),
              plotOutput("pathway_cnetplot", height = "760px")
            ),
            div(
              class = "panel pathway-results-panel",
              h2("Pathway Results"),
              tableOutput("pathway_results_table"),
              div(
                class = "download-grid",
                downloadButton("download_pathway_results", "Pathway Results CSV"),
                downloadButton("download_significant_pathways", "Significant Pathways CSV"),
                downloadButton("download_pathway_top_plot", "Top Plotted Pathways CSV"),
                downloadButton("download_pathway_ranks", "Ranked Genes CSV")
              )
            ),
            div(
              class = "panel",
              h2("Pathway Detail"),
              uiOutput("pathway_term_selector"),
              plotOutput("pathway_enrichment_plot", height = "360px"),
              h3("Leading-Edge Genes"),
              tableOutput("pathway_leading_edge_table"),
              downloadButton("download_pathway_leading_edge", "Leading-Edge Genes CSV")
            ),
            div(
              class = "panel",
              h2("Leading-Edge Expression"),
              uiOutput("pathway_heatmap_status"),
              plotOutput("pathway_heatmap", height = "540px")
            )
          ),
          tabPanel(
            "Network Analysis",
            div(
              class = "panel",
              h2("Weighted Gene Co-Expression Network Analysis"),
              uiOutput("wgcna_status"),
              div(
                class = "download-grid",
                numericInput("wgcna_max_genes", "Maximum variable genes", value = 5000, min = 10, max = 50000, step = 100),
                numericInput("wgcna_soft_power", "Soft-threshold power", value = 6, min = 1, max = 30, step = 1),
                numericInput("wgcna_min_module_size", "Minimum module size", value = 20, min = 2, max = 5000, step = 1)
              ),
              div(
                class = "download-grid",
                numericInput("wgcna_merge_cut_height", "Module merge cut height", value = 0.25, min = 0, max = 1, step = 0.05),
                selectInput(
                  "wgcna_network_type",
                  "Network type",
                  choices = c("Signed" = "signed", "Unsigned" = "unsigned"),
                  selected = "signed",
                  selectize = FALSE
                ),
                selectInput(
                  "wgcna_correlation",
                  "Correlation",
                  choices = c("Pearson" = "pearson", "Biweight midcorrelation" = "bicor"),
                  selected = "pearson",
                  selectize = FALSE
                )
              ),
              div(
                class = "download-grid",
                selectInput(
                  "wgcna_expression_scale",
                  "Normalized-expression scale",
                  choices = c(
                    "Use workflow scale" = "auto",
                    "Linear normalized values" = "linear",
                    "Already log2-transformed" = "log2"
                  ),
                  selected = "auto",
                  selectize = FALSE
                )
              ),
              actionButton(
                "run_wgcna",
                "Run WGCNA",
                class = "btn-primary run-button"
              ),
              tableOutput("wgcna_summary")
            ),
            div(
              class = "panel",
              h2("Module Sizes"),
              plotOutput("wgcna_module_plot", height = "420px")
            ),
            div(
              class = "panel",
              h2("Module-Trait Correlations"),
              uiOutput("wgcna_trait_status"),
              plotOutput("wgcna_trait_heatmap", height = "480px"),
              tableOutput("wgcna_trait_table")
            ),
            div(
              class = "panel",
              h2("Gene Modules"),
              tableOutput("wgcna_gene_table"),
              div(
                class = "download-grid",
                downloadButton("download_wgcna_module_summary", "Module Summary CSV"),
                downloadButton("download_wgcna_gene_modules", "Gene Modules CSV"),
                downloadButton("download_wgcna_traits", "Module-Trait CSV"),
                downloadButton("download_wgcna_eigengenes", "Module Eigengenes CSV")
              )
            )
          ),
          tabPanel(
            "Plots",
            div(
              class = "panel",
              h2("PCA"),
              plotOutput("pca_plot", height = "360px")
            ),
            div(
              class = "panel",
              h2("Volcano / Fold Change"),
              plotOutput("volcano_plot", height = "390px")
            ),
            div(
              class = "panel",
              h2("Regulation Summary"),
              plotOutput("regulation_summary_plot", height = "330px")
            ),
            div(
              class = "panel",
              h2("MA / Mean Expression"),
              plotOutput("ma_plot", height = "390px")
            )
          ),
          tabPanel(
            "Downloads",
            div(
              class = "panel",
              h2("Downloads"),
              h3("Current Results"),
              uiOutput("download_controls"),
              h3("Example Inputs"),
              uiOutput("example_download_controls")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  app_dir <- normalizePath(".", mustWork = TRUE)

  output$runtime_status <- renderText({
    deseq_version <- if (requireNamespace("DESeq2", quietly = TRUE)) {
      as.character(utils::packageVersion("DESeq2"))
    } else {
      "not installed"
    }
    fgsea_version <- if (requireNamespace("fgsea", quietly = TRUE)) {
      as.character(utils::packageVersion("fgsea"))
    } else {
      "not installed"
    }
    wgcna_version <- if (requireNamespace("WGCNA", quietly = TRUE)) {
      as.character(utils::packageVersion("WGCNA"))
    } else {
      "not installed"
    }
    sprintf("R %s | DESeq2 %s | fgsea %s | WGCNA %s", as.character(getRversion()), deseq_version, fgsea_version, wgcna_version)
  })

  output$metadata_upload <- renderUI({
    req(input$input_type)
    if (!input$input_type %in% matrix_modes) {
      return(div(class = "soft-note", "Fold-change modes do not need a metadata file."))
    }

    tagList(
      fileInput(
        "metadata_file",
        "Sample metadata (optional)",
        accept = c(".csv", ".tsv", ".txt")
      ),
      div(class = "soft-note", "If metadata is not supplied, groups are inferred from sample names such as condition_Rep1.")
    )
  })

  output$analysis_controls <- renderUI({
    req(input$input_type)

    matrix_controls <- if (input$input_type %in% matrix_modes) {
      tagList(
        selectInput("condition_col", "Condition column", choices = character(), selectize = FALSE),
        selectInput("reference_level", "Reference group", choices = character(), selectize = FALSE),
        selectInput("treatment_level", "Comparison group", choices = character(), selectize = FALSE)
      )
    }

    mode_controls <- switch(
      input$input_type,
      counts = tagList(
        selectInput("batch_col", "Primary adjustment/batch column", choices = c("None" = "__none__"), selectize = FALSE),
        selectizeInput(
          "adjustment_cols",
          "Additional factors to adjust for",
          choices = character(),
          multiple = TRUE,
          options = list(plugins = list("remove_button"))
        ),
        checkboxInput("use_interaction_design", "Use interaction design / custom contrast", value = FALSE),
        conditionalPanel(
          condition = "input.use_interaction_design == true",
          selectInput("interaction_col", "Interaction factor", choices = c("None" = "__none__"), selectize = FALSE),
          selectInput("interaction_reference_level", "Interaction reference level", choices = character(), selectize = FALSE),
          selectInput("interaction_comparison_level", "Interaction comparison level", choices = character(), selectize = FALSE),
          selectInput(
            "deseq_contrast_mode",
            "DESeq2 contrast to report",
            choices = c(
              "Main condition effect, adjusted for all factors" = "condition",
              "Interaction: difference in condition effect between factor levels" = "interaction",
              "Condition effect within selected interaction level" = "condition_at_interaction",
              "Custom DESeq2 result name" = "custom_name"
            ),
            selected = "interaction",
            selectize = FALSE
          ),
          textInput("custom_results_name", "Custom DESeq2 result name", value = ""),
          div(class = "soft-note", "Adjustment factors are included as additive terms. Interaction mode uses condition + interaction factor + condition:interaction factor. Custom names must match DESeq2 resultsNames().")
        ),
        checkboxInput("use_preprocessed_counts", "Use CPM-filtered counts for DESeq2", value = TRUE),
        numericInput("min_total_count", "Minimum total count after preprocessing", value = 0, min = 0, step = 1)
      ),
      normalized = tagList(
        selectInput(
          "expression_scale",
          "Expression scale",
          choices = c(
            "Auto-detect" = "auto",
            "Linear normalized values" = "linear",
            "Already log2-transformed" = "log2"
          ),
          selected = "auto",
          selectize = FALSE
        )
      ),
      fc_padj = tagList(
        selectInput(
          "fold_change_scale",
          "Fold-change scale",
          choices = c("Auto-detect" = "auto", "Linear fold change" = "linear", "Log2 fold change" = "log2"),
          selected = "auto",
          selectize = FALSE
        )
      ),
      fc_only = tagList(
        selectInput(
          "fold_change_scale",
          "Fold-change scale",
          choices = c("Auto-detect" = "auto", "Linear fold change" = "linear", "Log2 fold change" = "log2"),
          selected = "auto",
          selectize = FALSE
        )
      )
    )

    tagList(
      matrix_controls,
      mode_controls,
      numericInput("alpha", "Adjusted p-value threshold", value = 0.1, min = 0.001, max = 1, step = 0.01),
      numericInput("min_fold_change", "Minimum fold-change", value = 2, min = 1, step = 0.25),
      actionButton("run_analysis", "Run Analysis", class = "btn-primary run-button")
    )
  })

  effective_lfc_cutoff <- reactive({
    req(input$min_fold_change)
    log2(max(input$min_fold_change, 1))
  })

  selected_adjustment_cols <- reactive({
    if (!identical(input$input_type, "counts")) {
      return(NULL)
    }
    columns <- c(
      if (!is.null(input$batch_col) && !identical(input$batch_col, "__none__")) input$batch_col else NULL,
      if (is.null(input$adjustment_cols)) NULL else input$adjustment_cols
    )
    columns <- normalize_design_columns(columns)
    if (!is.null(input$condition_col)) {
      columns <- setdiff(columns, input$condition_col)
    }
    if (isTRUE(input$use_interaction_design) && !is.null(input$interaction_col) && !identical(input$interaction_col, "__none__")) {
      columns <- setdiff(columns, input$interaction_col)
    }
    columns
  })

  selected_interaction_col <- reactive({
    if (!identical(input$input_type, "counts") || !isTRUE(input$use_interaction_design)) {
      return(NULL)
    }
    if (is.null(input$interaction_col) || identical(input$interaction_col, "__none__")) {
      return(NULL)
    }
    input$interaction_col
  })

  selected_deseq_contrast_mode <- reactive({
    if (!identical(input$input_type, "counts") || !isTRUE(input$use_interaction_design)) {
      return("condition")
    }
    if (is.null(input$deseq_contrast_mode)) "condition" else input$deseq_contrast_mode
  })

  read_primary_raw <- reactive({
    req(input$input_type)
    if (isTRUE(input$use_example)) {
      example_path <- file.path(app_dir, "sample_data", example_file_for_mode(input$input_type))
      return(read_dge_table(example_path, basename(example_path)))
    }
    req(input$count_file)
    read_dge_table(input$count_file$datapath, input$count_file$name)
  })

  observe({
    req(identical(input$input_type, "counts"))
    raw_data <- read_primary_raw()
    count_check <- tryCatch(
      {
        prepare_count_matrix(raw_data)
        NULL
      },
      error = function(err) conditionMessage(err)
    )

    if (!is.null(count_check) && grepl("decimal values", count_check, fixed = TRUE)) {
      showNotification(
        "Decimal values detected. Switching to Normalized expression data mode.",
        type = "warning",
        duration = 8
      )
      updateSelectInput(session, "input_type", selected = "normalized")
    }
  })

  read_metadata_optional <- reactive({
    req(input$input_type)
    if (!input$input_type %in% matrix_modes) {
      return(NULL)
    }

    if (isTRUE(input$use_example)) {
      return(prepare_sample_metadata(read_dge_table(file.path(app_dir, "sample_data", "metadata.csv"), "metadata.csv")))
    }

    if (is.null(input$metadata_file)) {
      return(NULL)
    }

    prepare_sample_metadata(read_dge_table(input$metadata_file$datapath, input$metadata_file$name))
  })

  matrix_values <- reactive({
    req(input$input_type)
    if (identical(input$input_type, "counts")) {
      return(prepare_count_matrix(read_primary_raw()))
    }
    if (identical(input$input_type, "normalized")) {
      return(prepare_expression_matrix(read_primary_raw()))
    }
    NULL
  })

  matrix_metadata <- reactive({
    req(input$input_type)
    req(input$input_type %in% matrix_modes)
    metadata_for_matrix(matrix_values(), read_metadata_optional())
  })

  fold_change_preview <- reactive({
    req(input$input_type)
    req(input$input_type %in% c("fc_padj", "fc_only"))
    standardize_fold_change_table(
      read_primary_raw(),
      require_adjusted_p = identical(input$input_type, "fc_padj"),
      fold_change_scale = if (is.null(input$fold_change_scale)) "auto" else input$fold_change_scale
    )
  })

  preprocessed_counts <- reactive({
    req(input$input_type)
    req(identical(input$input_type, "counts"))
    preprocess_count_matrix(
      count_matrix = matrix_values(),
      metadata = matrix_metadata(),
      min_cpm = input$preprocess_min_cpm,
      min_samples = input$preprocess_min_samples,
      pseudocount = input$preprocess_pseudocount
    )
  })

  selected_annotation <- reactive({
    req(input$annotation_species)
    if (identical(input$annotation_species, "none")) {
      return(NULL)
    }
    read_gene_annotation(file.path(app_dir, "annotations"), input$annotation_species)
  })

  annotated_results <- reactive({
    result <- analysis_result()
    req(result)
    annotate_result_table(
      result_table = result$results,
      annotation_info = selected_annotation(),
      match_mode = if (is.null(input$annotation_match)) "auto" else input$annotation_match
    )
  })

  observe({
    req(input$input_type)
    req(input$input_type %in% matrix_modes)
    metadata <- matrix_metadata()
    choices <- analysis_columns(metadata)
    selected_condition <- choices[grepl("^condition$|^group$|^treatment$", tolower(choices))][1]
    if (is.na(selected_condition)) {
      selected_condition <- choices[1]
    }

    updateSelectInput(
      session,
      "condition_col",
      choices = choices,
      selected = selected_condition
    )

    if (identical(input$input_type, "counts")) {
      factor_choices <- setdiff(choices, selected_condition)
      updateSelectInput(
        session,
        "batch_col",
        choices = c("None" = "__none__", factor_choices),
        selected = "__none__"
      )
      updateSelectizeInput(
        session,
        "adjustment_cols",
        choices = factor_choices,
        selected = character(),
        server = FALSE
      )
      updateSelectInput(
        session,
        "interaction_col",
        choices = c("None" = "__none__", factor_choices),
        selected = "__none__"
      )
    }
  })

  observe({
    req(input$input_type)
    req(input$input_type %in% matrix_modes)
    metadata <- matrix_metadata()
    req(input$condition_col)
    req(input$condition_col %in% names(metadata))
    levels_found <- condition_levels(metadata, input$condition_col)
    selected_reference <- levels_found[1]
    selected_treatment <- if (length(levels_found) >= 2) levels_found[2] else levels_found[1]

    updateSelectInput(session, "reference_level", choices = levels_found, selected = selected_reference)
    updateSelectInput(session, "treatment_level", choices = levels_found, selected = selected_treatment)
  })

  observe({
    req(input$input_type)
    req(identical(input$input_type, "counts"))
    metadata <- matrix_metadata()
    req(input$condition_col)
    choices <- setdiff(analysis_columns(metadata), input$condition_col)

    current_batch <- if (!is.null(input$batch_col) && input$batch_col %in% choices) input$batch_col else "__none__"
    current_adjustments <- intersect(if (is.null(input$adjustment_cols)) character() else input$adjustment_cols, choices)
    current_interaction <- if (!is.null(input$interaction_col) && input$interaction_col %in% choices) input$interaction_col else "__none__"

    updateSelectInput(session, "batch_col", choices = c("None" = "__none__", choices), selected = current_batch)
    updateSelectizeInput(session, "adjustment_cols", choices = choices, selected = current_adjustments, server = FALSE)
    updateSelectInput(session, "interaction_col", choices = c("None" = "__none__", choices), selected = current_interaction)
  })

  observe({
    req(input$input_type)
    req(identical(input$input_type, "counts"))
    req(input$interaction_col)
    if (identical(input$interaction_col, "__none__")) {
      updateSelectInput(session, "interaction_reference_level", choices = character(), selected = character())
      updateSelectInput(session, "interaction_comparison_level", choices = character(), selected = character())
      return()
    }

    metadata <- matrix_metadata()
    req(input$interaction_col %in% names(metadata))
    levels_found <- condition_levels(metadata, input$interaction_col)
    selected_reference <- levels_found[1]
    selected_comparison <- if (length(levels_found) >= 2) levels_found[2] else levels_found[1]
    updateSelectInput(session, "interaction_reference_level", choices = levels_found, selected = selected_reference)
    updateSelectInput(session, "interaction_comparison_level", choices = levels_found, selected = selected_comparison)
  })

  output$input_status <- renderUI({
    loaded <- FALSE
    message <- "No input loaded."

    tryCatch(
      {
        if (input$input_type %in% matrix_modes) {
          values <- matrix_values()
          metadata <- matrix_metadata()
          loaded <- TRUE
          message <- sprintf(
            "%s genes, %s samples, %s metadata rows.",
            nrow(values),
            ncol(values),
            nrow(metadata)
          )
        } else {
          fc <- fold_change_preview()
          loaded <- TRUE
          message <- sprintf(
            "%s genes loaded. Fold-change column: %s.",
            nrow(fc$result_table),
            fc$fold_change_column
          )
        }
      },
      error = function(err) {
        message <<- conditionMessage(err)
      }
    )

    if (!loaded) {
      return(diagnostic_alert(
        message,
        severity = "warning",
        input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
        input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
        input_type = input$input_type,
        context = "input upload"
      ))
    }

    div(
      class = "alert-block",
      message
    )
  })

  output$annotation_status <- renderUI({
    if (is.null(input$annotation_species) || identical(input$annotation_species, "none")) {
      return(div(class = "soft-note", "No built-in annotation selected."))
    }

    tryCatch(
      {
        annotation <- selected_annotation()
        manifest <- annotation$manifest
        result <- analysis_result()
        if (is.null(result)) {
          return(div(
            class = "alert-block",
            sprintf(
              "%s loaded: %s genes, %s, %s.",
              manifest$label[1],
              format(nrow(annotation$table), big.mark = ","),
              manifest$dataset[1],
              manifest$assembly[1]
            )
          ))
        }

        annotated <- annotated_results()
        summary <- attr(annotated, "annotation_summary")
        div(
          class = "alert-block",
          sprintf(
            "%s matched %s of %s result genes by %s.",
            summary$label,
            format(summary$matched, big.mark = ","),
            format(summary$total, big.mark = ","),
            if (identical(summary$match_mode, "symbol")) "gene symbol" else "Ensembl/SGD ID"
          )
        )
      },
      error = function(err) {
        diagnostic_alert(
          conditionMessage(err),
          severity = "warning",
          input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
          input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
          input_type = input$input_type,
          context = "annotation"
        )
      }
    )
  })

  output$custom_gene_set_upload <- renderUI({
    if (!identical(input$gene_set_source, "custom")) {
      return(div(class = "soft-note", "Using the selected built-in GO, TF.Target.GTRD, or optional KEGG source. Custom GMT/CSV upload is optional. KEGG pathway mappings are not bundled; selecting KEGG downloads them from KEGG REST into a local user cache, so internet access and compliance with KEGG terms are required."))
    }

    tagList(
      fileInput(
        "geneset_file",
        "Custom gene set file",
        accept = c(".gmt", ".csv", ".tsv", ".txt")
      ),
      div(class = "soft-note", "Use this only for your own GMT or term-to-gene table. Do not upload the count matrix here.")
    )
  })

  observeEvent(input$enrichment_mode, {
    if (identical(input$enrichment_mode, "rosbys_lab")) {
      updateNumericInput(session, "ora_padj_cutoff", value = 0.01)
      updateNumericInput(session, "ora_top_n", value = 10)
    } else {
      updateNumericInput(session, "ora_padj_cutoff", value = 0.1)
      updateNumericInput(session, "ora_top_n", value = 20)
    }
  }, ignoreInit = TRUE)

  output$enrichment_mode_note <- renderUI({
    if (!identical(input$enrichment_mode, "rosbys_lab")) {
      return(NULL)
    }

    div(
      class = "soft-note",
      "Rosby's Lab-style ORA uses protein-coding genes as the background, tests upregulated and downregulated genes separately, keeps terms with at least 2 overlaps, applies FDR < 0.01 by default, and displays up to 5 terms per direction. It uses the selected GO, TF.Target.GTRD, KEGG, or custom gene-set source."
    )
  })

  output$preflight_metrics <- renderUI({
    tryCatch(
      {
        if (input$input_type %in% matrix_modes) {
          values <- matrix_values()
          metadata <- matrix_metadata()

          return(
            div(
              class = "metric-grid",
              div(class = "metric", span("Genes"), strong(format(nrow(values), big.mark = ","))),
              div(class = "metric", span("Samples"), strong(format(ncol(values), big.mark = ","))),
              div(class = "metric", span("Metadata rows"), strong(format(nrow(metadata), big.mark = ","))),
              div(class = "metric", span("Mode"), strong(if (input$input_type == "counts") "DESeq2" else "Expression"))
            )
          )
        }

        fc <- fold_change_preview()
        has_padj <- any(!is.na(fc$result_table$padj))
        div(
          class = "metric-grid",
          div(class = "metric", span("Genes"), strong(format(nrow(fc$result_table), big.mark = ","))),
          div(class = "metric", span("Fold-change scale"), strong(fc$fold_change_scale)),
          div(class = "metric", span("Adjusted p-values"), strong(if (has_padj) "Yes" else "No")),
          div(class = "metric", span("Mode"), strong("Table"))
        )
      },
      error = function(err) {
        diagnostic_alert(
          conditionMessage(err),
          severity = "warning",
          input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
          input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
          input_type = input$input_type,
          context = "preflight"
        )
      }
    )
  })

  output$preflight_table <- renderTable({
    tryCatch(
      {
        if (input$input_type %in% matrix_modes) {
          metadata <- matrix_metadata()
          req(input$condition_col)
          req(input$condition_col %in% names(metadata))

          condition_counts <- as.data.frame(table(metadata[[input$condition_col]]), stringsAsFactors = FALSE)
          names(condition_counts) <- c("group", "samples")
          return(condition_counts)
        }

        fc <- fold_change_preview()
        data.frame(
          field = c("fold_change_column", "adjusted_p_column"),
          value = c(fc$fold_change_column, if (is.null(fc$adjusted_p_column)) "not supplied" else fc$adjusted_p_column),
          stringsAsFactors = FALSE
        )
      },
      error = function(err) {
        issue <- diagnose_uploaded_table(
          input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
          input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
          input_type = input$input_type
        )
        if (is.null(issue) || !nzchar(issue)) {
          issue <- conditionMessage(err)
        }
        if (is.null(issue) || !nzchar(issue)) {
          issue <- "Input could not be summarized. Check the diagnostic card above."
        }
        data.frame(field = "Input issue", value = issue, stringsAsFactors = FALSE)
      }
    )
  })

  output$preflight_warnings <- renderUI({
    tryCatch(
      {
        warnings <- character()
        info_notes <- character()

        if (input$input_type %in% matrix_modes) {
          values <- matrix_values()
          metadata <- matrix_metadata()
          req(input$condition_col, input$treatment_level, input$reference_level)

          batch_col <- if (identical(input$batch_col, "__none__") || !identical(input$input_type, "counts")) NULL else input$batch_col
          adjustment_cols <- if (identical(input$input_type, "counts")) selected_adjustment_cols() else NULL
          interaction_col <- if (identical(input$input_type, "counts")) selected_interaction_col() else NULL
          warnings <- attr(metadata, "dge_warnings")
          validation_warnings <- tryCatch(
            validate_analysis_settings(
              abs(round(values)),
              metadata,
              input$condition_col,
              input$treatment_level,
              input$reference_level,
              batch_col,
              adjustment_cols = adjustment_cols,
              interaction_col = interaction_col,
              min_total_count = if (identical(input$input_type, "counts")) input$min_total_count else 0
            ),
            error = function(err) conditionMessage(err)
          )
          warnings <- unique(c(warnings, validation_warnings))

          if (identical(input$input_type, "counts")) {
            design_note <- sprintf(
              "DESeq2 design: %s.",
              if (!is.null(interaction_col)) {
                sprintf(
                  "~ %scondition + interaction + condition:interaction",
                  if (length(adjustment_cols) > 0) "adjustment factors + " else ""
                )
              } else if (length(adjustment_cols) > 0) {
                "~ adjustment factors + condition"
              } else {
                "~ condition"
              }
            )
            info_notes <- c(info_notes, design_note)
            if (!is.null(interaction_col)) {
              info_notes <- c(
                info_notes,
                "Interaction contrasts test whether the condition effect changes across the selected interaction factor; additive adjustment factors control for those effects in the DESeq2 model."
              )
            } else if (length(adjustment_cols) > 0) {
              info_notes <- c(info_notes, "Selected adjustment factors are included as additive DESeq2 model terms to control for those effects.")
            }
          }

          if (identical(input$input_type, "normalized")) {
            info_notes <- c(
              info_notes,
              "Ready: decimal expression values are allowed in Normalized expression mode. This mode computes group means, log2 fold changes, Welch t-test p-values, and BH-adjusted p-values."
            )
          }
        } else {
          fc <- fold_change_preview()
          warnings <- c(
            "Fold-change mode accepts result tables and does not run DESeq2.",
            if (!any(!is.na(fc$result_table$padj))) "No adjusted p-values are available, so significance calls are disabled."
          )
        }

        if (length(warnings) == 0) {
          return(div(
            class = "alert-block",
            tags$div("Inputs passed the current checks."),
            if (length(info_notes) > 0) tags$ul(lapply(unique(info_notes), tags$li))
          ))
        }

        div(
          class = "alert-block warning",
          tags$ul(lapply(warnings, tags$li)),
          if (length(info_notes) > 0) tags$ul(lapply(unique(info_notes), tags$li))
        )
      },
      error = function(err) {
        diagnostic_alert(
          conditionMessage(err),
          severity = "warning",
          input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
          input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
          input_type = input$input_type,
          context = "preflight"
        )
      }
    )
  })

  output$preprocess_status <- renderUI({
    if (!identical(input$input_type, "counts")) {
      return(div(
        class = "alert-block warning",
        "Preprocessing is available for Read counts data. Normalized expression and fold-change inputs are already downstream data types."
      ))
    }

    tryCatch(
      {
        processed <- preprocessed_counts()
        div(
          class = "alert-block",
          sprintf(
            "Ready: %s of %s genes pass CPM filtering. Filtered raw counts can be used for DESeq2; log2(CPM + c) is for QC/exploration.",
            format(processed$genes_after, big.mark = ","),
            format(processed$genes_before, big.mark = ",")
          )
        )
      },
      error = function(err) {
        diagnostic_alert(
          conditionMessage(err),
          severity = "warning",
          input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
          input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
          input_type = input$input_type,
          context = "preprocessing"
        )
      }
    )
  })

  output$preprocess_metrics <- renderUI({
    req(identical(input$input_type, "counts"))
    processed <- preprocessed_counts()
    div(
      class = "metric-grid",
      div(class = "metric", span("Genes before"), strong(format(processed$genes_before, big.mark = ","))),
      div(class = "metric", span("Genes after"), strong(format(processed$genes_after, big.mark = ","))),
      div(class = "metric", span("Genes removed"), strong(format(processed$genes_removed, big.mark = ","))),
      div(class = "metric", span("Filter"), strong(sprintf("%s CPM", processed$min_cpm)))
    )
  })

  output$preprocess_sample_table <- renderTable({
    req(identical(input$input_type, "counts"))
    processed <- preprocessed_counts()
    summary_table <- processed$sample_summary
    numeric_cols <- vapply(summary_table, is.numeric, logical(1))
    summary_table[numeric_cols] <- lapply(summary_table[numeric_cols], function(column) {
      format(round(column), big.mark = ",", scientific = FALSE)
    })
    summary_table
  })

  analysis_result <- reactiveVal(NULL)
  analysis_error <- reactiveVal(NULL)
  pathway_result <- reactiveVal(NULL)
  pathway_error <- reactiveVal(NULL)
  wgcna_result <- reactiveVal(NULL)
  wgcna_error <- reactiveVal(NULL)

  reset_app_state <- function() {
    analysis_result(NULL)
    analysis_error(NULL)
    pathway_result(NULL)
    pathway_error(NULL)
    wgcna_result(NULL)
    wgcna_error(NULL)
    updateTabsetPanel(session, "main_tabs", selected = "Preflight")
  }

  observeEvent(input$refresh_app, {
    reset_app_state()
    session$sendCustomMessage("transcriptoscope-refresh", list(reset = TRUE))
  }, ignoreInit = TRUE)

  significant_pathway_results <- function(pathway) {
    req(pathway)
    if (is.null(pathway$results) || !"padj" %in% names(pathway$results)) {
      return(pathway$results[FALSE, , drop = FALSE])
    }
    cutoff <- pathway$padj_cutoff
    if (is.null(cutoff) || length(cutoff) != 1 || is.na(cutoff)) {
      cutoff <- if (is.null(input$pathway_padj_cutoff)) 0.05 else input$pathway_padj_cutoff
    }
    pathway$results[!is.na(pathway$results$padj) & pathway$results$padj < cutoff, , drop = FALSE]
  }

  top_plotted_pathway_results <- function(pathway) {
    req(pathway)
    top_n <- if (is.null(input$pathway_top_n) || is.na(input$pathway_top_n)) 20 else input$pathway_top_n
    select_top_pathway_results(pathway$results, top_n = top_n, padj_cutoff = pathway$padj_cutoff)
  }

  go_ontology <- local({
    cache <- NULL
    function() {
      if (is.null(cache)) {
        cache <<- read_go_ontology(file.path(app_dir, "gene_sets"))
      }
      cache
    }
  })

  go_dag_top_n <- function() {
    if (is.null(input$go_dag_top_n) || is.na(input$go_dag_top_n)) {
      return(15)
    }
    input$go_dag_top_n
  }

  go_dag_ancestor_depth <- function() {
    if (is.null(input$go_dag_ancestor_depth) || is.na(input$go_dag_ancestor_depth)) {
      return(4)
    }
    input$go_dag_ancestor_depth
  }

  go_dag_padj_cutoff <- function() {
    if (is.null(input$go_dag_padj_cutoff) || is.na(input$go_dag_padj_cutoff)) {
      return(0.1)
    }
    input$go_dag_padj_cutoff
  }

  go_dag_data_for <- function(ora) {
    req(ora)
    build_go_dag(
      ora_table = ora$results,
      go_ontology = go_ontology(),
      top_n = go_dag_top_n(),
      padj_cutoff = go_dag_padj_cutoff(),
      max_ancestor_depth = go_dag_ancestor_depth()
    )
  }

  pathway_cnet_top_n <- function() {
    if (is.null(input$pathway_cnet_top_n) || is.na(input$pathway_cnet_top_n)) {
      return(5)
    }
    input$pathway_cnet_top_n
  }

  pathway_cnet_max_genes <- function() {
    if (is.null(input$pathway_cnet_max_genes) || is.na(input$pathway_cnet_max_genes)) {
      return(15)
    }
    input$pathway_cnet_max_genes
  }

  pathway_cnet_edges_for <- function(pathway) {
    req(pathway)
    pathway_cnetplot_edges(
      pathway,
      top_n = pathway_cnet_top_n(),
      max_genes_per_pathway = pathway_cnet_max_genes(),
      padj_cutoff = pathway$padj_cutoff,
      show_ids = isTRUE(input$pathway_show_ids)
    )
  }

  save_pathway_cnetplot_file <- function(file, pathway, format = c("png", "tiff"), dpi = 300) {
    format <- match.arg(format)
    plot <- make_pathway_cnetplot(
      pathway,
      top_n = pathway_cnet_top_n(),
      max_genes_per_pathway = pathway_cnet_max_genes(),
      padj_cutoff = pathway$padj_cutoff,
      show_ids = isTRUE(input$pathway_show_ids)
    )
    args <- list(
      filename = file,
      plot = plot,
      width = 12,
      height = 8.5,
      units = "in",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
    if (identical(format, "png")) {
      args$device <- "png"
    } else {
      args$device <- "tiff"
      args$compression <- "lzw"
    }
    do.call(ggplot2::ggsave, args)
  }

  observeEvent(input$run_analysis, {
    analysis_result(NULL)
    analysis_error(NULL)
    pathway_result(NULL)
    pathway_error(NULL)
    wgcna_result(NULL)
    wgcna_error(NULL)

    withProgress(message = "Running analysis", value = 0.1, {
      tryCatch(
        {
          mode <- input$input_type
          incProgress(0.25, detail = "Checking inputs")

          result <- switch(
            mode,
            counts = run_deseq2_workflow(
              count_matrix = if (isTRUE(input$use_preprocessed_counts)) {
                preprocessed_counts()$filtered_counts
              } else {
                matrix_values()
              },
              metadata = matrix_metadata(),
              condition_col = input$condition_col,
              treatment_level = input$treatment_level,
              reference_level = input$reference_level,
              batch_col = if (identical(input$batch_col, "__none__")) NULL else input$batch_col,
              adjustment_cols = if (is.null(input$adjustment_cols)) NULL else input$adjustment_cols,
              interaction_col = selected_interaction_col(),
              interaction_reference_level = if (is.null(input$interaction_reference_level)) NULL else input$interaction_reference_level,
              interaction_comparison_level = if (is.null(input$interaction_comparison_level)) NULL else input$interaction_comparison_level,
              contrast_mode = selected_deseq_contrast_mode(),
              custom_results_name = if (is.null(input$custom_results_name)) NULL else input$custom_results_name,
              min_total_count = input$min_total_count,
              alpha = input$alpha
            ),
            normalized = run_normalized_expression_workflow(
              expression_matrix = matrix_values(),
              metadata = read_metadata_optional(),
              condition_col = input$condition_col,
              treatment_level = input$treatment_level,
              reference_level = input$reference_level,
              expression_scale = input$expression_scale,
              alpha = input$alpha
            ),
            fc_padj = run_fold_change_workflow(
              data = read_primary_raw(),
              require_adjusted_p = TRUE,
              fold_change_scale = input$fold_change_scale,
              alpha = input$alpha
            ),
            fc_only = run_fold_change_workflow(
              data = read_primary_raw(),
              require_adjusted_p = FALSE,
              fold_change_scale = input$fold_change_scale,
              alpha = input$alpha
            )
          )

          incProgress(0.65, detail = "Preparing results")
          analysis_result(result)
          updateTabsetPanel(session, "main_tabs", selected = "Results")
        },
        error = function(err) {
          analysis_error(conditionMessage(err))
        }
      )
    })
  })

  output$analysis_status <- renderUI({
    err <- analysis_error()
    if (!is.null(err)) {
      return(diagnostic_alert(
        err,
        severity = "error",
        input_name = if (!is.null(input$count_file)) input$count_file$name else NULL,
        input_path = if (!is.null(input$count_file)) input$count_file$datapath else NULL,
        input_type = input$input_type,
        context = "analysis"
      ))
    }

    result <- analysis_result()
    if (is.null(result)) {
      ready_message <- tryCatch(
        {
          if (input$input_type %in% matrix_modes) {
            values <- matrix_values()
            metadata <- matrix_metadata()
            sprintf("Ready to run: %s genes and %s samples loaded. Click Run Analysis.", nrow(values), ncol(values))
          } else {
            fc <- fold_change_preview()
            sprintf("Ready to run: %s gene rows loaded. Click Run Analysis.", nrow(fc$result_table))
          }
        },
        error = function(err) NULL
      )

      if (!is.null(ready_message)) {
        return(div(class = "alert-block", ready_message))
      }

      return(div(class = "alert-block warning", "Load data, check the preflight tab, then run the analysis."))
    }

    has_padj <- any(!is.na(result$results$padj))
    if (has_padj) {
      sig <- sum(!is.na(result$results$padj) & result$results$padj < result$alpha)
      message <- sprintf(
        "%s complete: %s genes processed, %s genes pass adjusted p < %s.",
        result$workflow_label,
        format(result$genes_tested, big.mark = ","),
        format(sig, big.mark = ","),
        result$alpha
      )
    } else {
      message <- sprintf(
        "%s complete: %s genes processed. No adjusted p-values were supplied.",
        result$workflow_label,
        format(result$genes_tested, big.mark = ",")
      )
    }

    div(class = "alert-block", message)
  })

  output$results_table <- renderTable({
    annotated <- annotated_results()
    req(annotated)
    format_result_preview(annotated, 100)
  })

  output$pca_plot <- renderPlot({
    result <- analysis_result()
    req(result)
    validate(need(!is.null(result$pca_data), "PCA is available for matrix inputs only."))
    make_pca_plot(result)
  })

  output$volcano_plot <- renderPlot({
    result <- analysis_result()
    req(result)
    make_volcano_plot(result$results, result$alpha, effective_lfc_cutoff())
  })

  output$regulation_summary_plot <- renderPlot({
    result <- analysis_result()
    req(result)
    make_regulation_summary_plot(result$results, result$alpha, effective_lfc_cutoff())
  })

  output$ma_plot <- renderPlot({
    result <- analysis_result()
    req(result)
    draw_ma_plot(result)
  })

  gene_sets_for_enrichment <- reactive({
    result <- analysis_result()
    result_table <- annotated_results()
    req(result, result_table)

    if (identical(input$gene_set_source, "custom")) {
      req(input$geneset_file)
      custom_sets <- read_gene_set_file(input$geneset_file$datapath, input$geneset_file$name)
      attr(custom_sets, "gene_set_summary") <- list(
        label = sprintf("Custom file: %s", input$geneset_file$name),
        dataset = "custom",
        assembly = "",
        domain = "custom",
        match_mode = "custom",
        matched = NA_integer_,
        total = nrow(result_table),
        terms = length(unique(custom_sets$term_id)),
        rows = nrow(custom_sets)
      )
      return(custom_sets)
    }

    collection_info <- gene_set_collection_info(input$go_domain, default = "go_all")
    collection <- collection_info$collection
    domain <- collection_info$domain

    selected_annotation_key <- input$gene_set_source
    if (is.null(selected_annotation_key) || identical(selected_annotation_key, "auto")) {
      req(input$annotation_species)
      selected_annotation_key <- input$annotation_species
    }
    matched <- gene_set_manifest[
      gene_set_manifest$annotation_key == selected_annotation_key &
        gene_set_manifest$collection == collection,
      ,
      drop = FALSE
    ]
    validate(need(nrow(matched) > 0, "Choose an organism with the selected gene set collection."))

    built_in <- read_builtin_gene_sets(
      gene_set_dir = file.path(app_dir, "gene_sets"),
      gene_set_key = matched$key[1],
      domain = domain
    )
    prepare_builtin_gene_sets_for_results(
      result_table = result_table,
      built_in_gene_sets = built_in,
      match_mode = if (is.null(input$annotation_match)) "auto" else input$annotation_match
    )
  })

  selected_enrichment_genes <- reactive({
    result <- analysis_result()
    req(result)
    select_de_genes(
      result_table = result$results,
      alpha = result$alpha,
      lfc_cutoff = effective_lfc_cutoff(),
      gene_list = if (is.null(input$ora_gene_list)) "both" else input$ora_gene_list
    )
  })

  enrichment_result <- reactive({
    result <- analysis_result()
    req(result)
    gene_sets <- gene_sets_for_enrichment()

    if (identical(input$enrichment_mode, "rosbys_lab")) {
      return(run_rosbys_lab_enrichment(
        result_table = result$results,
        gene_sets = gene_sets,
        annotation_info = selected_annotation(),
        alpha = result$alpha,
        lfc_cutoff = effective_lfc_cutoff(),
        padj_cutoff = if (is.null(input$ora_padj_cutoff)) 0.01 else input$ora_padj_cutoff,
        max_terms_per_direction = 5,
        min_overlap = 2,
        case_sensitive = identical(input$gene_set_source, "custom") && isTRUE(input$ora_case_sensitive)
      ))
    }

    selected <- selected_enrichment_genes()

    run_ora_analysis(
      universe_genes = result$results$gene_id,
      selected_genes = selected$genes,
      gene_sets = gene_sets,
      min_set_size = if (is.null(input$ora_min_set_size)) 3 else input$ora_min_set_size,
      max_set_size = if (is.null(input$ora_max_set_size)) 2000 else input$ora_max_set_size,
      padj_cutoff = if (is.null(input$ora_padj_cutoff)) 0.1 else input$ora_padj_cutoff,
      case_sensitive = identical(input$gene_set_source, "custom") && isTRUE(input$ora_case_sensitive)
    )
  })

  output$enrichment_status <- renderUI({
    result <- analysis_result()
    if (is.null(result)) {
      return(div(class = "alert-block warning", "Run the DGE analysis first. The enrichment tab uses the current result table automatically."))
    }
    if (identical(input$gene_set_source, "custom") && is.null(input$geneset_file)) {
      return(div(class = "alert-block warning", "Upload a custom GMT/CSV/TSV gene set file, or choose a built-in GO, TF.Target.GTRD, or optional KEGG source."))
    }

    tryCatch(
      {
        gene_sets <- gene_sets_for_enrichment()
        gene_set_summary <- attr(gene_sets, "gene_set_summary")
        ora <- enrichment_result()
        database_note <- if (!is.null(gene_set_summary)) {
          sprintf(
            " Database: %s; matched %s of %s result genes by %s.",
            gene_set_summary$label,
            if (is.na(gene_set_summary$matched)) "custom" else format(gene_set_summary$matched, big.mark = ","),
            format(gene_set_summary$total, big.mark = ","),
            if (identical(gene_set_summary$match_mode, "symbol")) "gene symbol" else if (identical(gene_set_summary$match_mode, "ensembl")) "Ensembl/SGD ID" else "custom IDs"
          )
        } else {
          ""
        }

        if (identical(ora$mode, "rosbys_lab")) {
          up_count <- sum(ora$selected_table$regulation == "Upregulated", na.rm = TRUE)
          down_count <- sum(ora$selected_table$regulation == "Downregulated", na.rm = TRUE)
          tested <- ora$summary$value[match("Gene sets tested", ora$summary$metric)]
          note <- sprintf(
            "Rosby's Lab-style ORA ready: %s upregulated and %s downregulated protein-coding DEG genes, %s tested gene sets, %s displayed terms at enrichment FDR < %s.%s",
            format(up_count, big.mark = ","),
            format(down_count, big.mark = ","),
            format(as.numeric(tested), big.mark = ","),
            format(nrow(ora$results), big.mark = ","),
            ora$padj_cutoff,
            database_note
          )
        } else {
          selected <- selected_enrichment_genes()
          significant_terms <- sum(!is.na(ora$results$padj) & ora$results$padj < ora$padj_cutoff)
          note <- if (selected$has_padj) {
            sprintf(
              "ORA ready: %s selected DEG genes, %s tested gene sets, %s terms pass FDR < %s.%s",
              format(selected$selected_count, big.mark = ","),
              format(nrow(ora$results), big.mark = ","),
              format(significant_terms, big.mark = ","),
              ora$padj_cutoff,
              database_note
            )
          } else {
            sprintf(
              "ORA ready: %s fold-change-selected genes, %s tested gene sets, %s terms pass FDR < %s.%s No adjusted p-values were available in the DGE input.",
              format(selected$selected_count, big.mark = ","),
              format(nrow(ora$results), big.mark = ","),
              format(significant_terms, big.mark = ","),
              ora$padj_cutoff,
              database_note
            )
          }
        }
        div(class = "alert-block", note)
      },
      error = function(err) {
        diagnostic_alert(
          conditionMessage(err),
          severity = "warning",
          input_name = if (!is.null(input$geneset_file)) input$geneset_file$name else NULL,
          input_type = input$input_type,
          context = "ORA enrichment"
        )
      }
    )
  })

  output$enrichment_summary <- renderTable({
    ora <- enrichment_result()
    ora$summary
  })

  output$enrichment_plot <- renderPlot({
    ora <- enrichment_result()
    make_ora_plot(
      ora$results,
      top_n = if (is.null(input$ora_top_n)) 20 else input$ora_top_n,
      padj_cutoff = if (is.null(input$ora_padj_cutoff)) 0.1 else input$ora_padj_cutoff
    )
  })

  output$go_dag_status <- renderUI({
    if (is.null(analysis_result())) {
      return(div(class = "alert-block warning", "Run the DGE analysis first. The GO DAG uses the current ORA enrichment result."))
    }
    ora <- tryCatch(enrichment_result(), error = function(err) err)
    if (inherits(ora, "error")) {
      return(diagnostic_alert(
        conditionMessage(ora),
        severity = "warning",
        input_name = if (!is.null(input$geneset_file)) input$geneset_file$name else NULL,
        input_type = input$input_type,
        context = "GO DAG"
      ))
    }
    go_rows <- ora$results[grepl("^GO:[0-9]{7}$", ora$results$term_id), , drop = FALSE]
    if (nrow(go_rows) == 0) {
      return(div(class = "alert-block warning", "The GO DAG plot is available for enrichment results with GO term IDs. Choose a GO collection, or upload custom sets whose term IDs are GO accessions."))
    }
    cutoff <- go_dag_padj_cutoff()
    significant <- sum(!is.na(go_rows$padj) & go_rows$padj < cutoff)
    if (significant == 0) {
      return(div(class = "alert-block warning", sprintf("No over-represented GO terms pass DAG FDR < %s.", cutoff)))
    }
    ontology_status <- tryCatch(
      {
        ontology <- go_ontology()
        sprintf(
          " Ontology: %s GO terms and %s relationships loaded.",
          format(nrow(ontology$terms), big.mark = ","),
          format(nrow(ontology$edges), big.mark = ",")
        )
      },
      error = function(err) conditionMessage(err)
    )
    div(
      class = "soft-note",
      sprintf(
        "DAG ready: %s over-represented GO terms pass FDR < %s; the graph will draw up to %s top terms plus %s ancestor level(s).%s",
        format(significant, big.mark = ","),
        cutoff,
        go_dag_top_n(),
        go_dag_ancestor_depth(),
        ontology_status
      )
    )
  })

  output$go_dag_plot <- renderPlot({
    ora <- enrichment_result()
    go_rows <- ora$results[grepl("^GO:[0-9]{7}$", ora$results$term_id), , drop = FALSE]
    if (nrow(go_rows) == 0) {
      graphics::plot.new()
      graphics::text(
        0.5,
        0.5,
        "Choose a GO collection to draw the GO DAG plot.",
        cex = 1.1
      )
      return(invisible(NULL))
    }
    significant <- go_rows[!is.na(go_rows$padj) & go_rows$padj < go_dag_padj_cutoff(), , drop = FALSE]
    if (nrow(significant) == 0) {
      graphics::plot.new()
      graphics::text(
        0.5,
        0.5,
        sprintf("No GO terms pass DAG FDR < %s.", go_dag_padj_cutoff()),
        cex = 1.1
      )
      return(invisible(NULL))
    }
    make_go_dag_plot(
      ora_table = ora$results,
      go_ontology = go_ontology(),
      top_n = go_dag_top_n(),
      padj_cutoff = go_dag_padj_cutoff(),
      max_ancestor_depth = go_dag_ancestor_depth()
    )
  })

  output$enrichment_table <- renderTable({
    ora <- enrichment_result()
    if (identical(ora$mode, "rosbys_lab")) {
      format_rosbys_lab_enrichment_preview(ora$results, 100)
    } else {
      format_enrichment_preview(ora$results, 100)
    }
  })

  pathway_gene_sets <- reactive({
    result_table <- annotated_results()
    req(result_table)

    collection_info <- gene_set_collection_info(input$pathway_collection, default = "kegg_pathway")
    collection <- collection_info$collection
    domain <- collection_info$domain

    selected_annotation_key <- input$pathway_gene_set_source
    if (is.null(selected_annotation_key) || identical(selected_annotation_key, "auto")) {
      req(input$annotation_species)
      selected_annotation_key <- input$annotation_species
    }
    validate(need(!identical(selected_annotation_key, "none"), "Choose an organism annotation or a pathway organism."))
    matched <- gene_set_manifest[
      gene_set_manifest$annotation_key == selected_annotation_key &
        gene_set_manifest$collection == collection,
      ,
      drop = FALSE
    ]
    validate(need(nrow(matched) > 0, "Choose an organism with the selected pathway collection."))

    built_in <- read_builtin_gene_sets(
      gene_set_dir = file.path(app_dir, "gene_sets"),
      gene_set_key = matched$key[1],
      domain = domain
    )
    prepare_builtin_gene_sets_for_results(
      result_table = result_table,
      built_in_gene_sets = built_in,
      match_mode = if (is.null(input$annotation_match)) "auto" else input$annotation_match
    )
  })

  observeEvent(input$run_pathway_analysis, {
    pathway_result(NULL)
    pathway_error(NULL)
    result <- analysis_result()
    if (is.null(result)) {
      pathway_error("Run the DGE analysis before running pathway analysis.")
      return()
    }

    withProgress(message = "Running ranked pathway analysis", value = 0.1, {
      tryCatch(
        {
          incProgress(0.25, detail = "Matching ranked genes to pathways")
          gene_sets <- pathway_gene_sets()
          result_table <- annotated_results()
          pathway <- run_preranked_pathway_analysis(
            result_table = result_table,
            gene_sets = gene_sets,
            rank_metric = if (is.null(input$pathway_rank_metric)) "log2fc" else input$pathway_rank_metric,
            min_set_size = if (is.null(input$pathway_min_set_size)) 3 else input$pathway_min_set_size,
            max_set_size = if (is.null(input$pathway_max_set_size)) 500 else input$pathway_max_set_size,
            padj_cutoff = if (is.null(input$pathway_padj_cutoff)) 0.05 else input$pathway_padj_cutoff,
            gene_padj_cutoff = if (is.null(input$pathway_gene_padj_cutoff)) 1 else input$pathway_gene_padj_cutoff,
            absolute_ranking = isTRUE(input$pathway_absolute_ranking)
          )
          attr(pathway, "gene_set_summary") <- attr(gene_sets, "gene_set_summary")
          incProgress(0.6, detail = "Preparing pathway plots and leading-edge genes")
          pathway_result(pathway)
        },
        error = function(err) {
          pathway_error(conditionMessage(err))
        }
      )
    })
  })

  output$pathway_status <- renderUI({
    if (is.null(analysis_result())) {
      return(div(class = "alert-block warning", "Run the DGE analysis first. This tab automatically ranks the current result genes."))
    }
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      return(diagnostic_alert(
        "The fgsea package is not installed. Run Install_Packages.bat, then restart TranscriptoScope.",
        severity = "error",
        input_type = input$input_type,
        context = "pathway analysis"
      ))
    }
    err <- pathway_error()
    if (!is.null(err)) {
      return(diagnostic_alert(
        err,
        severity = "warning",
        input_type = input$input_type,
        context = "pathway analysis"
      ))
    }
    pathway <- pathway_result()
    if (is.null(pathway)) {
      return(div(
        class = "alert-block",
        "Ready for preranked GSEA. The full DGE result is used by default; set Maximum gene FDR included below 1 only when you deliberately want to filter the ranking."
      ))
    }

    gene_set_summary <- attr(pathway, "gene_set_summary")
    significant <- sum(!is.na(pathway$results$padj) & pathway$results$padj < pathway$padj_cutoff)
    database_note <- if (!is.null(gene_set_summary)) {
      sprintf(
        " Database: %s; matched %s of %s result genes by %s.",
        gene_set_summary$label,
        format(gene_set_summary$matched, big.mark = ","),
        format(gene_set_summary$total, big.mark = ","),
        if (identical(gene_set_summary$match_mode, "symbol")) "gene symbol" else "Ensembl/SGD ID"
      )
    } else {
      ""
    }
    div(
      class = "alert-block",
      sprintf(
        "Pathway analysis complete: %s ranked genes, %s pathways tested, and %s pathways pass FDR < %s.%s",
        format(nrow(pathway$ranked_table), big.mark = ","),
        format(nrow(pathway$results), big.mark = ","),
        format(significant, big.mark = ","),
        pathway$padj_cutoff,
        database_note
      )
    )
  })

  output$pathway_summary <- renderTable({
    pathway <- pathway_result()
    req(pathway)
    pathway$summary
  })

  output$pathway_summary_plot <- renderPlot({
    pathway <- pathway_result()
    req(pathway)
    make_pathway_summary_plot(
      pathway$results,
      top_n = if (is.null(input$pathway_top_n)) 20 else input$pathway_top_n,
      padj_cutoff = pathway$padj_cutoff,
      show_ids = isTRUE(input$pathway_show_ids)
    )
  })

  output$pathway_cnetplot <- renderPlot({
    pathway <- pathway_result()
    req(pathway)
    make_pathway_cnetplot(
      pathway,
      top_n = pathway_cnet_top_n(),
      max_genes_per_pathway = pathway_cnet_max_genes(),
      padj_cutoff = pathway$padj_cutoff,
      show_ids = isTRUE(input$pathway_show_ids)
    )
  })

  output$pathway_results_table <- renderTable({
    pathway <- pathway_result()
    req(pathway)
    format_pathway_preview(pathway$results, 100, show_ids = isTRUE(input$pathway_show_ids))
  })

  output$pathway_term_selector <- renderUI({
    pathway <- pathway_result()
    req(pathway)
    choices <- stats::setNames(
      pathway$results$term_id,
      sprintf(
        "%s | NES %.2f | FDR %s",
        pathway$results$term_name,
        pathway$results$NES,
        formatC(pathway$results$padj, format = "e", digits = 2)
      )
    )
    selectInput(
      "pathway_term",
      "Pathway",
      choices = choices,
      selected = pathway$results$term_id[1],
      selectize = TRUE
    )
  })

  output$pathway_enrichment_plot <- renderPlot({
    pathway <- pathway_result()
    req(pathway, input$pathway_term)
    req(input$pathway_term %in% pathway$results$term_id)
    make_pathway_enrichment_plot(pathway, input$pathway_term)
  })

  output$pathway_leading_edge_table <- renderTable({
    pathway <- pathway_result()
    req(pathway, input$pathway_term)
    req(input$pathway_term %in% pathway$results$term_id)
    leading <- pathway_leading_edge_table(pathway, input$pathway_term)
    leading$pvalue <- ifelse(is.na(leading$pvalue), "NA", formatC(leading$pvalue, format = "e", digits = 3))
    leading$padj <- ifelse(is.na(leading$padj), "NA", formatC(leading$padj, format = "e", digits = 3))
    leading$rank_score <- round(leading$rank_score, 4)
    leading$log2FoldChange <- round(leading$log2FoldChange, 4)
    leading
  })

  output$pathway_heatmap_status <- renderUI({
    result <- analysis_result()
    pathway <- pathway_result()
    if (is.null(pathway)) {
      return(div(class = "soft-note", "Run pathway analysis and choose a pathway to display its leading-edge expression."))
    }
    if (is.null(result$normalized_counts)) {
      return(div(class = "alert-block warning", "A heatmap needs sample-level counts or expression data. Fold-change-only inputs still provide pathway scores and leading-edge genes."))
    }
    div(class = "soft-note", "Rows are standardized within each gene. Read-count workflows use log2(normalized count + 1) before row scaling.")
  })

  output$pathway_heatmap <- renderPlot({
    result <- analysis_result()
    pathway <- pathway_result()
    req(result, pathway, input$pathway_term)
    req(input$pathway_term %in% pathway$results$term_id)
    validate(need(!is.null(result$normalized_counts), "Sample-level expression data are not available for this input mode."))
    leading <- pathway_leading_edge_table(pathway, input$pathway_term)
    make_pathway_heatmap(result, leading$gene_id, max_genes = 40)
  })

  observeEvent(input$run_wgcna, {
    wgcna_result(NULL)
    wgcna_error(NULL)
    result <- analysis_result()
    if (is.null(result)) {
      wgcna_error("Run the DGE analysis before running WGCNA.")
      return()
    }
    if (is.null(result$normalized_counts)) {
      wgcna_error("WGCNA needs sample-level count or expression data. Fold-change-only inputs cannot be used.")
      return()
    }
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      wgcna_error("The WGCNA package is not installed. Run Install_Packages.bat, then restart TranscriptoScope.")
      return()
    }

    withProgress(message = "Running WGCNA", value = 0.1, {
      tryCatch(
        {
          incProgress(0.25, detail = "Filtering variable genes")
          result <- run_wgcna_analysis(
            workflow_result = result,
            max_genes = if (is.null(input$wgcna_max_genes)) 5000 else input$wgcna_max_genes,
            min_module_size = if (is.null(input$wgcna_min_module_size)) 20 else input$wgcna_min_module_size,
            soft_power = if (is.null(input$wgcna_soft_power)) 6 else input$wgcna_soft_power,
            merge_cut_height = if (is.null(input$wgcna_merge_cut_height)) 0.25 else input$wgcna_merge_cut_height,
            network_type = if (is.null(input$wgcna_network_type)) "signed" else input$wgcna_network_type,
            cor_type = if (is.null(input$wgcna_correlation)) "pearson" else input$wgcna_correlation,
            expression_scale = if (is.null(input$wgcna_expression_scale)) "auto" else input$wgcna_expression_scale
          )
          incProgress(0.65, detail = "Preparing module tables")
          wgcna_result(result)
        },
        error = function(err) {
          wgcna_error(conditionMessage(err))
        }
      )
    })
  })

  output$wgcna_status <- renderUI({
    result <- analysis_result()
    if (is.null(result)) {
      return(div(class = "alert-block warning", "Run the DGE analysis first. WGCNA uses the current sample-level matrix."))
    }
    if (is.null(result$normalized_counts)) {
      return(div(class = "alert-block warning", "WGCNA needs sample-level count or expression data. Fold-change-only inputs cannot be used."))
    }
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      return(diagnostic_alert(
        "The WGCNA package is not installed. Run Install_Packages.bat, then restart TranscriptoScope.",
        severity = "error",
        input_type = input$input_type,
        context = "WGCNA"
      ))
    }
    err <- wgcna_error()
    if (!is.null(err)) {
      return(diagnostic_alert(
        err,
        severity = "warning",
        input_type = input$input_type,
        context = "WGCNA"
      ))
    }
    wgcna <- wgcna_result()
    if (is.null(wgcna)) {
      return(div(class = "alert-block", sprintf(
        "Ready for WGCNA: %s genes and %s samples are available.",
        format(result$genes_tested, big.mark = ","),
        format(result$samples_tested, big.mark = ",")
      )))
    }

    module_count <- sum(wgcna$module_summary$module != "grey")
    div(
      class = "alert-block",
      sprintf(
        "WGCNA complete: %s genes, %s samples, and %s assigned module(s).",
        format(wgcna$settings$genes_used, big.mark = ","),
        format(wgcna$settings$samples_used, big.mark = ","),
        format(module_count, big.mark = ",")
      )
    )
  })

  output$wgcna_summary <- renderTable({
    wgcna <- wgcna_result()
    req(wgcna)
    format_wgcna_module_summary(wgcna$module_summary, 100)
  })

  output$wgcna_module_plot <- renderPlot({
    wgcna <- wgcna_result()
    req(wgcna)
    make_wgcna_module_plot(wgcna)
  })

  output$wgcna_trait_status <- renderUI({
    wgcna <- wgcna_result()
    if (is.null(wgcna)) {
      return(div(class = "soft-note", "Run WGCNA to calculate module-trait correlations."))
    }
    if (nrow(wgcna$trait_correlations) == 0) {
      return(div(class = "alert-block warning", "No usable metadata traits were available for module-trait correlations."))
    }
    div(class = "soft-note", "Module eigengenes are correlated with numeric metadata and encoded categorical metadata.")
  })

  output$wgcna_trait_heatmap <- renderPlot({
    wgcna <- wgcna_result()
    req(wgcna)
    validate(need(nrow(wgcna$trait_correlations) > 0, "No module-trait correlations are available to plot."))
    make_wgcna_trait_heatmap(wgcna)
  })

  output$wgcna_trait_table <- renderTable({
    wgcna <- wgcna_result()
    req(wgcna)
    validate(need(nrow(wgcna$trait_correlations) > 0, "No module-trait correlations are available."))
    format_wgcna_trait_correlations(wgcna$trait_correlations, 100)
  })

  output$wgcna_gene_table <- renderTable({
    wgcna <- wgcna_result()
    req(wgcna)
    format_wgcna_gene_modules(wgcna$gene_modules, 100)
  })

  output$download_controls <- renderUI({
    controls <- list(
      downloadButton("download_results", "Results CSV"),
      downloadButton("download_matrix", if (input$input_type %in% matrix_modes) "Matrix Values" else "Processed Table"),
      downloadButton("download_bundle", "Result Bundle")
    )
    if (!is.null(input$annotation_species) && !identical(input$annotation_species, "none")) {
      controls <- c(controls, list(downloadButton("download_annotation", "Annotation CSV")))
    }
    do.call(div, c(list(class = "download-grid"), controls))
  })

  output$example_download_controls <- renderUI({
    controls <- lapply(example_download_specs, function(spec) downloadButton(spec$id, spec$label))
    do.call(div, c(list(class = "download-grid"), controls))
  })

  copy_example_input <- function(example_name, destination) {
    source <- file.path(app_dir, "sample_data", example_name)
    validate(need(file.exists(source), sprintf("The bundled example file %s is missing.", example_name)))
    file.copy(source, destination, overwrite = TRUE)
  }

  for (spec in example_download_specs) {
    local({
      spec <- spec
      output[[spec$id]] <- downloadHandler(
        filename = function() {
          if (identical(spec$id, "download_example_bundle")) {
            sprintf("transcriptoscope_example_inputs_%s.zip", format(Sys.Date(), "%Y%m%d"))
          } else {
            sprintf("transcriptoscope_%s", spec$file)
          }
        },
        content = function(file) {
          if (!identical(spec$id, "download_example_bundle")) {
            copy_example_input(spec$file, file)
            return(invisible(TRUE))
          }

          example_dir <- file.path(app_dir, "sample_data")
          validate(need(dir.exists(example_dir), "The bundled sample_data folder is missing."))
          bundle_dir <- file.path(tempdir(), sprintf("transcriptoscope_examples_%s", as.integer(Sys.time())))
          dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
          example_files <- stats::na.omit(vapply(example_download_specs, `[[`, character(1), "file"))
          missing_examples <- example_files[!file.exists(file.path(example_dir, example_files))]
          validate(need(length(missing_examples) == 0, paste("Missing bundled example file(s):", paste(missing_examples, collapse = ", "))))
          file.copy(file.path(example_dir, example_files), bundle_dir, overwrite = TRUE)
          writeLines(
            c(
              "TranscriptoScope example input files",
              "",
              "counts.csv: raw integer read-count matrix for DESeq2 mode.",
              "metadata.csv: sample metadata matching the count and expression matrix column names.",
              "normalized_expression.csv: decimal expression matrix for normalized-expression mode.",
              "fold_changes_padj.csv: table for fold changes and adjusted p-values mode.",
              "fold_changes_only.csv: table for fold-changes-only mode.",
              "gene_sets.csv: custom term-to-gene example for ORA or ranked pathway uploads.",
              "",
              "Keep the first row as headers. Keep gene IDs in the first column of matrix-style files."
            ),
            file.path(bundle_dir, "README_example_inputs.txt")
          )
          oldwd <- setwd(bundle_dir)
          on.exit(setwd(oldwd), add = TRUE)
          utils::zip(zipfile = file, files = list.files(bundle_dir), flags = "-q")
          invisible(TRUE)
        }
      )
    })
  }

  output$download_results <- downloadHandler(
    filename = function() {
      sprintf("dge_results_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      annotated <- annotated_results()
      req(annotated)
      utils::write.csv(annotated, file, row.names = FALSE)
    }
  )

  output$download_annotation <- downloadHandler(
    filename = function() {
      annotation <- selected_annotation()
      req(annotation)
      sprintf("ensembl_annotation_%s_%s.csv", annotation$manifest$key[1], format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      annotation <- selected_annotation()
      req(annotation)
      utils::write.csv(annotation$table, file, row.names = FALSE)
    }
  )

  output$download_enrichment <- downloadHandler(
    filename = function() {
      if (identical(input$enrichment_mode, "rosbys_lab")) {
        sprintf("rosbys_lab_style_enrichment_%s.csv", format(Sys.Date(), "%Y%m%d"))
      } else {
        sprintf("ora_enrichment_%s.csv", format(Sys.Date(), "%Y%m%d"))
      }
    },
    content = function(file) {
      ora <- enrichment_result()
      req(ora)
      if (identical(ora$mode, "rosbys_lab")) {
        utils::write.csv(format_rosbys_lab_enrichment_export(ora$results), file, row.names = FALSE)
      } else {
        utils::write.csv(ora$results, file, row.names = FALSE)
      }
    }
  )

  output$download_selected_ora_genes <- downloadHandler(
    filename = function() {
      if (identical(input$enrichment_mode, "rosbys_lab")) {
        sprintf("rosbys_lab_style_selected_genes_%s.csv", format(Sys.Date(), "%Y%m%d"))
      } else {
        sprintf("ora_selected_genes_%s.csv", format(Sys.Date(), "%Y%m%d"))
      }
    },
    content = function(file) {
      ora <- tryCatch(enrichment_result(), error = function(err) NULL)
      if (!is.null(ora) && identical(ora$mode, "rosbys_lab")) {
        utils::write.csv(ora$selected_table, file, row.names = FALSE)
      } else {
        selected <- selected_enrichment_genes()
        req(selected)
        utils::write.csv(selected$table, file, row.names = FALSE)
      }
    }
  )

  output$download_go_dag_plot <- downloadHandler(
    filename = function() {
      sprintf("go_dag_plot_%s.png", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      ora <- enrichment_result()
      req(ora)
      ggplot2::ggsave(
        filename = file,
        plot = make_go_dag_plot(
          ora_table = ora$results,
          go_ontology = go_ontology(),
          top_n = go_dag_top_n(),
          padj_cutoff = go_dag_padj_cutoff(),
          max_ancestor_depth = go_dag_ancestor_depth()
        ),
        width = 14,
        height = 8,
        dpi = 300,
        bg = "white"
      )
    }
  )

  output$download_go_dag_nodes <- downloadHandler(
    filename = function() {
      sprintf("go_dag_nodes_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      ora <- enrichment_result()
      req(ora)
      utils::write.csv(go_dag_data_for(ora)$nodes, file, row.names = FALSE)
    }
  )

  output$download_go_dag_edges <- downloadHandler(
    filename = function() {
      sprintf("go_dag_edges_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      ora <- enrichment_result()
      req(ora)
      utils::write.csv(go_dag_data_for(ora)$edges, file, row.names = FALSE)
    }
  )

  output$download_pathway_results <- downloadHandler(
    filename = function() {
      sprintf("pathway_analysis_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      utils::write.csv(pathway$results, file, row.names = FALSE)
    }
  )

  output$download_significant_pathways <- downloadHandler(
    filename = function() {
      sprintf("significant_pathways_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      utils::write.csv(significant_pathway_results(pathway), file, row.names = FALSE)
    }
  )

  output$download_pathway_top_plot <- downloadHandler(
    filename = function() {
      sprintf("top_plotted_pathways_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      utils::write.csv(top_plotted_pathway_results(pathway), file, row.names = FALSE)
    }
  )

  output$download_pathway_cnetplot <- downloadHandler(
    filename = function() {
      sprintf("pathway_cnetplot_%s.png", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      save_pathway_cnetplot_file(file, pathway, format = "png", dpi = 300)
    }
  )

  output$download_pathway_cnetplot_tiff <- downloadHandler(
    filename = function() {
      sprintf("pathway_cnetplot_%s.tiff", format(Sys.Date(), "%Y%m%d"))
    },
    contentType = "image/tiff",
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      save_pathway_cnetplot_file(file, pathway, format = "tiff", dpi = 600)
    }
  )

  output$download_pathway_cnet_edges <- downloadHandler(
    filename = function() {
      sprintf("pathway_cnetplot_edges_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      utils::write.csv(pathway_cnet_edges_for(pathway), file, row.names = FALSE)
    }
  )

  output$download_pathway_ranks <- downloadHandler(
    filename = function() {
      sprintf("pathway_ranked_genes_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway)
      utils::write.csv(pathway$ranked_table, file, row.names = FALSE)
    }
  )

  output$download_pathway_leading_edge <- downloadHandler(
    filename = function() {
      term <- if (is.null(input$pathway_term)) "selected_pathway" else gsub("[^A-Za-z0-9_-]", "_", input$pathway_term)
      sprintf("leading_edge_%s_%s.csv", term, format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      pathway <- pathway_result()
      req(pathway, input$pathway_term)
      req(input$pathway_term %in% pathway$results$term_id)
      utils::write.csv(pathway_leading_edge_table(pathway, input$pathway_term), file, row.names = FALSE)
    }
  )

  output$download_wgcna_module_summary <- downloadHandler(
    filename = function() {
      sprintf("wgcna_module_summary_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      wgcna <- wgcna_result()
      req(wgcna)
      utils::write.csv(wgcna$module_summary, file, row.names = FALSE)
    }
  )

  output$download_wgcna_gene_modules <- downloadHandler(
    filename = function() {
      sprintf("wgcna_gene_modules_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      wgcna <- wgcna_result()
      req(wgcna)
      utils::write.csv(wgcna$gene_modules, file, row.names = FALSE)
    }
  )

  output$download_wgcna_traits <- downloadHandler(
    filename = function() {
      sprintf("wgcna_module_trait_correlations_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      wgcna <- wgcna_result()
      req(wgcna)
      utils::write.csv(wgcna$trait_correlations, file, row.names = FALSE)
    }
  )

  output$download_wgcna_eigengenes <- downloadHandler(
    filename = function() {
      sprintf("wgcna_module_eigengenes_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      wgcna <- wgcna_result()
      req(wgcna)
      utils::write.csv(wgcna_module_eigengene_table(wgcna), file, row.names = FALSE)
    }
  )

  output$download_preprocessed_counts <- downloadHandler(
    filename = function() {
      sprintf("filtered_counts_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      processed <- preprocessed_counts()
      utils::write.csv(matrix_to_export_table(processed$filtered_counts), file, row.names = FALSE)
    }
  )

  output$download_preprocessed_expression <- downloadHandler(
    filename = function() {
      sprintf("log2_cpm_expression_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      processed <- preprocessed_counts()
      utils::write.csv(matrix_to_export_table(processed$transformed_expression), file, row.names = FALSE)
    }
  )

  output$download_preprocessed_metadata <- downloadHandler(
    filename = function() {
      sprintf("sample_metadata_%s.csv", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      processed <- preprocessed_counts()
      utils::write.csv(processed$metadata, file, row.names = FALSE)
    }
  )

  output$download_matrix <- downloadHandler(
    filename = function() {
      if (input$input_type %in% matrix_modes) {
        sprintf("matrix_values_%s.csv", format(Sys.Date(), "%Y%m%d"))
      } else {
        sprintf("processed_table_%s.csv", format(Sys.Date(), "%Y%m%d"))
      }
    },
    content = function(file) {
      result <- analysis_result()
      req(result)
      export_table <- if (is.null(result$normalized_counts)) annotated_results() else result$normalized_counts
      utils::write.csv(export_table, file, row.names = FALSE)
    }
  )

  output$download_bundle <- downloadHandler(
    filename = function() {
      sprintf("dge_result_bundle_%s.zip", format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      result <- analysis_result()
      req(result)
      annotated <- annotated_results()
      uploaded_name <- function(upload) {
        if (is.null(upload) || is.null(upload$name) || !nzchar(upload$name)) {
          return(NA_character_)
        }
        upload$name
      }

      bundle_dir <- file.path(tempdir(), sprintf("dge_bundle_%s", as.integer(Sys.time())))
      dir.create(bundle_dir, recursive = TRUE)
      export_result <- result
      export_result$results <- annotated

      annotation <- tryCatch(selected_annotation(), error = function(err) NULL)
      ora <- tryCatch(enrichment_result(), error = function(err) NULL)
      selected <- tryCatch(selected_enrichment_genes(), error = function(err) NULL)
      pathway <- pathway_result()
      wgcna <- wgcna_result()
      enrichment_gene_sets_used <- tryCatch(gene_sets_for_enrichment(), error = function(err) NULL)
      enrichment_gene_set_summary <- if (is.null(enrichment_gene_sets_used)) NULL else attr(enrichment_gene_sets_used, "gene_set_summary")
      pathway_gene_sets_used <- if (is.null(pathway)) NULL else tryCatch(pathway_gene_sets(), error = function(err) NULL)
      pathway_gene_set_summary <- if (is.null(pathway_gene_sets_used)) {
        if (is.null(pathway)) NULL else attr(pathway, "gene_set_summary")
      } else {
        attr(pathway_gene_sets_used, "gene_set_summary")
      }
      app_version <- tryCatch(readLines(file.path(app_dir, "VERSION"), warn = FALSE, n = 1), error = function(err) "unknown")
      bundle_reproducibility <- list(
        generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        app_version = app_version,
        input = list(
          input_type = input$input_type,
          input_file = if (isTRUE(input$use_example)) example_file_for_mode(input$input_type) else uploaded_name(input$count_file),
          metadata_file = if (isTRUE(input$use_example)) "bundled example or inferred metadata" else uploaded_name(input$metadata_file),
          used_bundled_example = isTRUE(input$use_example)
        ),
        analysis = list(
          condition_column = input$condition_col,
          reference_group = input$reference_level,
          comparison_group = input$treatment_level,
          batch_column = if (identical(input$batch_col, "__none__")) "none" else input$batch_col,
          adjustment_columns = paste(selected_adjustment_cols(), collapse = ";"),
          interaction_column = if (is.null(selected_interaction_col())) "none" else selected_interaction_col(),
          interaction_reference_level = if (is.null(input$interaction_reference_level)) NA_character_ else input$interaction_reference_level,
          interaction_comparison_level = if (is.null(input$interaction_comparison_level)) NA_character_ else input$interaction_comparison_level,
          deseq_contrast_mode = selected_deseq_contrast_mode(),
          custom_deseq_result_name = if (is.null(input$custom_results_name)) NA_character_ else input$custom_results_name,
          fitted_design_formula = if (is.null(export_result$design_formula)) NA_character_ else export_result$design_formula,
          deseq_result_contrast = if (is.null(export_result$result_contrast_used)) NA_character_ else export_result$result_contrast_used,
          deseq_result_name = if (is.null(export_result$result_name_used)) NA_character_ else export_result$result_name_used,
          adjusted_p_value_threshold = input$alpha,
          minimum_fold_change = input$min_fold_change,
          absolute_log2_fold_change_cutoff = effective_lfc_cutoff(),
          expression_scale = if (is.null(input$expression_scale)) NA_character_ else input$expression_scale,
          fold_change_scale = if (is.null(input$fold_change_scale)) NA_character_ else input$fold_change_scale
        ),
        preprocessing = list(
          minimum_cpm = input$preprocess_min_cpm,
          minimum_samples = input$preprocess_min_samples,
          transform_pseudocount = input$preprocess_pseudocount,
          use_cpm_filtered_counts_for_deseq2 = isTRUE(input$use_preprocessed_counts),
          minimum_total_count_after_preprocessing = if (identical(input$input_type, "counts")) input$min_total_count else NA_real_
        ),
        annotation = list(
          species = input$annotation_species,
          match_mode = input$annotation_match,
          annotation_rows_exported = if (is.null(annotation)) 0 else nrow(annotation$table)
        ),
        enrichment = if (is.null(ora)) {
          list(status = "not run")
        } else {
          list(
            status = "run",
            method = input$enrichment_mode,
            gene_set_organism = input$gene_set_source,
            gene_set_collection = input$go_domain,
            gene_set_file = uploaded_name(input$geneset_file),
            selected_deg_list = input$ora_gene_list,
            minimum_gene_set_size = input$ora_min_set_size,
            maximum_gene_set_size = input$ora_max_set_size,
            fdr_cutoff = input$ora_padj_cutoff,
            top_terms_to_plot = input$ora_top_n,
            case_sensitive_custom_matching = isTRUE(input$ora_case_sensitive),
            gene_set_label = if (is.null(enrichment_gene_set_summary)) NA_character_ else enrichment_gene_set_summary$label
          )
        },
        pathway = if (is.null(pathway)) {
          list(status = "not run")
        } else {
          list(
            status = "run",
            gene_set_organism = input$pathway_gene_set_source,
            pathway_collection = input$pathway_collection,
            ranking_metric = input$pathway_rank_metric,
            minimum_pathway_size = input$pathway_min_set_size,
            maximum_pathway_size = input$pathway_max_set_size,
            pathway_fdr_cutoff = input$pathway_padj_cutoff,
            maximum_gene_fdr_included = input$pathway_gene_padj_cutoff,
            use_absolute_ranking_scores = isTRUE(input$pathway_absolute_ranking),
            top_pathways_to_plot = input$pathway_top_n,
            show_pathway_ids = isTRUE(input$pathway_show_ids),
            gene_set_label = if (is.null(pathway_gene_set_summary)) NA_character_ else pathway_gene_set_summary$label
          )
        },
        wgcna = if (is.null(wgcna)) {
          list(status = "not run")
        } else {
          c(
            list(status = "run"),
            wgcna$settings
          )
        }
      )

      helper_source <- file.path(app_dir, "R", "deseq_helpers.R")
      if (file.exists(helper_source)) {
        file.copy(helper_source, file.path(bundle_dir, "transcriptoscope_deseq_helpers.R"), overwrite = TRUE)
      }
      write_result_bundle(export_result, bundle_dir, effective_lfc_cutoff(), reproducibility = bundle_reproducibility)

      if (!is.null(annotation)) {
        utils::write.csv(annotation$table, file.path(bundle_dir, "ensembl_annotation.csv"), row.names = FALSE)
      }

      if (!is.null(ora)) {
        if (identical(ora$mode, "rosbys_lab")) {
          utils::write.csv(format_rosbys_lab_enrichment_export(ora$results), file.path(bundle_dir, "rosbys_lab_style_enrichment.csv"), row.names = FALSE)
          utils::write.csv(ora$selected_table, file.path(bundle_dir, "rosbys_lab_style_selected_genes.csv"), row.names = FALSE)
        } else {
          utils::write.csv(ora$results, file.path(bundle_dir, "ora_enrichment.csv"), row.names = FALSE)
        }
        go_dag <- tryCatch(go_dag_data_for(ora), error = function(err) NULL)
        if (!is.null(go_dag) && nrow(go_dag$nodes) > 0) {
          utils::write.csv(go_dag$nodes, file.path(bundle_dir, "go_dag_nodes.csv"), row.names = FALSE)
          utils::write.csv(go_dag$edges, file.path(bundle_dir, "go_dag_edges.csv"), row.names = FALSE)
          tryCatch(
            ggplot2::ggsave(
              filename = file.path(bundle_dir, "go_dag_plot.png"),
              plot = make_go_dag_plot(
                ora_table = ora$results,
                go_ontology = go_ontology(),
                top_n = go_dag_top_n(),
                padj_cutoff = go_dag_padj_cutoff(),
                max_ancestor_depth = go_dag_ancestor_depth()
              ),
              width = 14,
              height = 8,
              dpi = 220,
              bg = "white"
            ),
            error = function(err) NULL
          )
        }
      }
      if (!is.null(selected) && (is.null(ora) || !identical(ora$mode, "rosbys_lab"))) {
        utils::write.csv(selected$table, file.path(bundle_dir, "ora_selected_genes.csv"), row.names = FALSE)
      }
      if (!is.null(enrichment_gene_sets_used)) {
        utils::write.csv(enrichment_gene_sets_used, file.path(bundle_dir, "enrichment_gene_sets_used.csv"), row.names = FALSE)
      }
      if (!is.null(ora) && identical(ora$mode, "rosbys_lab") && !is.null(ora$protein_coding_summary)) {
        utils::write.csv(
          data.frame(gene_id = ora$protein_coding_summary$universe_ids, stringsAsFactors = FALSE),
          file.path(bundle_dir, "rosbys_lab_protein_coding_universe.csv"),
          row.names = FALSE
        )
      }
      if (!is.null(pathway)) {
        utils::write.csv(pathway$results, file.path(bundle_dir, "pathway_analysis.csv"), row.names = FALSE)
        utils::write.csv(significant_pathway_results(pathway), file.path(bundle_dir, "pathway_significant_pathways.csv"), row.names = FALSE)
        utils::write.csv(top_plotted_pathway_results(pathway), file.path(bundle_dir, "pathway_top_plotted_pathways.csv"), row.names = FALSE)
        utils::write.csv(pathway$ranked_table, file.path(bundle_dir, "pathway_ranked_genes.csv"), row.names = FALSE)
        cnet_edges <- tryCatch(pathway_cnet_edges_for(pathway), error = function(err) NULL)
        if (!is.null(cnet_edges) && nrow(cnet_edges) > 0) {
          utils::write.csv(cnet_edges, file.path(bundle_dir, "pathway_cnetplot_edges.csv"), row.names = FALSE)
          tryCatch(
            save_pathway_cnetplot_file(file.path(bundle_dir, "pathway_cnetplot.png"), pathway, format = "png", dpi = 300),
            error = function(err) NULL
          )
          tryCatch(
            save_pathway_cnetplot_file(file.path(bundle_dir, "pathway_cnetplot.tiff"), pathway, format = "tiff", dpi = 600),
            error = function(err) NULL
          )
        }
      }
      if (!is.null(pathway_gene_sets_used)) {
        utils::write.csv(pathway_gene_sets_used, file.path(bundle_dir, "pathway_gene_sets_used.csv"), row.names = FALSE)
      }
      if (!is.null(wgcna)) {
        utils::write.csv(wgcna$module_summary, file.path(bundle_dir, "wgcna_module_summary.csv"), row.names = FALSE)
        utils::write.csv(wgcna$gene_modules, file.path(bundle_dir, "wgcna_gene_modules.csv"), row.names = FALSE)
        utils::write.csv(wgcna$trait_correlations, file.path(bundle_dir, "wgcna_module_trait_correlations.csv"), row.names = FALSE)
        utils::write.csv(wgcna_module_eigengene_table(wgcna), file.path(bundle_dir, "wgcna_module_eigengenes.csv"), row.names = FALSE)
      }
      write_analysis_report(export_result, bundle_dir, effective_lfc_cutoff(), bundle_reproducibility)

      oldwd <- setwd(bundle_dir)
      on.exit(setwd(oldwd), add = TRUE)
      utils::zip(zipfile = file, files = list.files(bundle_dir), flags = "-q")
    }
  )
}

shinyApp(ui, server)
