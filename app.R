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
  error = function(err) "0.4.4"
)
if (!nzchar(app_version) || is.na(app_version)) {
  app_version <- "0.4.4"
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
  "KEGG pathways" = "kegg_pathway"
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
        window.location.reload();
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
        p("Windows-friendly transcriptomics analysis for raw counts, expression matrices, fold-change tables, enrichment, and pathways"),
        p(class = "creator-line", "Created by Dr. Abubakar Abdulkadir | Dr. Rosby's Lab, Southern University A and M")
      ),
      div(
        class = "header-actions",
        actionButton("refresh_app", "Refresh", class = "refresh-button", title = "Reload the app"),
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
              uiOutput("download_controls")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  app_dir <- normalizePath(".", mustWork = TRUE)

  observeEvent(input$refresh_app, {
    session$sendCustomMessage("transcriptoscope-refresh", list())
  }, ignoreInit = TRUE)

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
    sprintf("R %s | DESeq2 %s | fgsea %s", as.character(getRversion()), deseq_version, fgsea_version)
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
        selectInput("batch_col", "Batch column", choices = c("None" = "__none__"), selectize = FALSE),
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
      updateSelectInput(
        session,
        "batch_col",
        choices = c("None" = "__none__", choices),
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

    div(
      class = if (loaded) "alert-block" else "alert-block warning",
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
        div(class = "alert-block warning", conditionMessage(err))
      }
    )
  })

  output$custom_gene_set_upload <- renderUI({
    if (!identical(input$gene_set_source, "custom")) {
      return(div(class = "soft-note", "Using the selected built-in GO or optional KEGG source. Custom GMT/CSV upload is optional. KEGG pathway mappings are not bundled; selecting KEGG downloads them from KEGG REST into a local user cache, so internet access and compliance with KEGG terms are required."))
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
      "Rosby's Lab-style ORA uses protein-coding genes as the background, tests upregulated and downregulated genes separately, keeps terms with at least 2 overlaps, applies FDR < 0.01 by default, and displays up to 5 terms per direction. It uses the selected GO, KEGG, or custom gene-set source."
    )
  })

  output$preflight_metrics <- renderUI({
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
  })

  output$preflight_table <- renderTable({
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
  })

  output$preflight_warnings <- renderUI({
    warnings <- character()
    info_notes <- character()

    if (input$input_type %in% matrix_modes) {
      values <- matrix_values()
      metadata <- matrix_metadata()
      req(input$condition_col, input$treatment_level, input$reference_level)

      batch_col <- if (identical(input$batch_col, "__none__") || !identical(input$input_type, "counts")) NULL else input$batch_col
      warnings <- attr(metadata, "dge_warnings")
      validation_warnings <- tryCatch(
        validate_analysis_settings(
          abs(round(values)),
          metadata,
          input$condition_col,
          input$treatment_level,
          input$reference_level,
          batch_col,
          if (identical(input$input_type, "counts")) input$min_total_count else 0
        ),
        error = function(err) conditionMessage(err)
      )
      warnings <- unique(c(warnings, validation_warnings))

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
        div(class = "alert-block warning", conditionMessage(err))
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

  observeEvent(input$run_analysis, {
    analysis_result(NULL)
    analysis_error(NULL)
    pathway_result(NULL)
    pathway_error(NULL)

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
      return(div(class = "alert-block error", err))
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
    req(result)

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
        total = nrow(result$results),
        terms = length(unique(custom_sets$term_id)),
        rows = nrow(custom_sets)
      )
      return(custom_sets)
    }

    collection_choice <- if (is.null(input$go_domain)) "go_all" else input$go_domain
    collection <- if (identical(collection_choice, "kegg_pathway")) "kegg" else "go"
    domain <- switch(
      collection_choice,
      go_biological_process = "biological_process",
      go_molecular_function = "molecular_function",
      go_cellular_component = "cellular_component",
      kegg_pathway = "kegg_pathway",
      "all"
    )

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
      result_table = result$results,
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
      return(div(class = "alert-block warning", "Upload a custom GMT/CSV/TSV gene set file, or choose a built-in GO collection or optional KEGG source."))
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
        div(class = "alert-block warning", conditionMessage(err))
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

  output$enrichment_table <- renderTable({
    ora <- enrichment_result()
    if (identical(ora$mode, "rosbys_lab")) {
      format_rosbys_lab_enrichment_preview(ora$results, 100)
    } else {
      format_enrichment_preview(ora$results, 100)
    }
  })

  pathway_gene_sets <- reactive({
    result <- analysis_result()
    req(result)

    collection_choice <- if (is.null(input$pathway_collection)) "kegg_pathway" else input$pathway_collection
    collection <- if (identical(collection_choice, "kegg_pathway")) "kegg" else "go"
    domain <- switch(
      collection_choice,
      go_biological_process = "biological_process",
      go_molecular_function = "molecular_function",
      go_cellular_component = "cellular_component",
      kegg_pathway = "kegg_pathway",
      "all"
    )

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
      result_table = result$results,
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
          pathway <- run_preranked_pathway_analysis(
            result_table = result$results,
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
      return(div(class = "alert-block error", "The fgsea package is not installed. Run Install_Packages.bat, then restart TranscriptoScope."))
    }
    err <- pathway_error()
    if (!is.null(err)) {
      return(div(class = "alert-block warning", err))
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
      }
      if (!is.null(pathway_gene_sets_used)) {
        utils::write.csv(pathway_gene_sets_used, file.path(bundle_dir, "pathway_gene_sets_used.csv"), row.names = FALSE)
      }
      write_analysis_report(export_result, bundle_dir, effective_lfc_cutoff(), bundle_reproducibility)

      oldwd <- setwd(bundle_dir)
      on.exit(setwd(oldwd), add = TRUE)
      utils::zip(zipfile = file, files = list.files(bundle_dir), flags = "-q")
    }
  )
}

shinyApp(ui, server)
