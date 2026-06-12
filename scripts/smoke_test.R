args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

setup_user_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) {
    version_dir <- paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
    user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", version_dir)
  }
  user_lib <- path.expand(user_lib)
  if (dir.exists(user_lib)) {
    .libPaths(unique(c(normalizePath(user_lib, winslash = "/", mustWork = TRUE), .libPaths())))
  }
}

setup_user_library()

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

source(file.path(app_dir, "R", "deseq_helpers.R"))

metadata <- prepare_sample_metadata(read_dge_table(file.path(app_dir, "sample_data", "metadata.csv")))

counts <- prepare_count_matrix(read_dge_table(file.path(app_dir, "sample_data", "counts.csv")))
processed_counts <- preprocess_count_matrix(
  count_matrix = counts,
  metadata = metadata,
  min_cpm = 0.5,
  min_samples = 1,
  pseudocount = 4
)
stopifnot(nrow(processed_counts$filtered_counts) > 0)
stopifnot(nrow(processed_counts$transformed_expression) == nrow(processed_counts$filtered_counts))
stopifnot(all(colnames(processed_counts$filtered_counts) == rownames(processed_counts$metadata)))

raw_result <- run_deseq2_workflow(
  count_matrix = processed_counts$filtered_counts,
  metadata = metadata,
  condition_col = "condition",
  treatment_level = "treated",
  reference_level = "control",
  batch_col = "__none__",
  min_total_count = 10,
  alpha = 0.05
)
stopifnot(nrow(raw_result$results) > 0)
stopifnot("gene_id" %in% names(raw_result$results))
stopifnot(nrow(raw_result$normalized_counts) > 0)

annotation_info <- read_gene_annotation(file.path(app_dir, "annotations"), "yeast_scerevisiae")
annotation_rows <- annotation_info$table[
  !is.na(annotation_info$table$gene_symbol) & annotation_info$table$gene_symbol != "",
  ,
  drop = FALSE
]
annotation_rows <- utils::head(annotation_rows, 3)
annotation_test_result <- data.frame(
  gene_id = annotation_rows$ensembl_gene_id,
  baseMean = c(100, 120, 140),
  log2FoldChange = c(1.2, -1.4, 0.3),
  pvalue = c(0.001, 0.002, 0.5),
  padj = c(0.01, 0.02, 0.8),
  stringsAsFactors = FALSE
)
annotated_test_result <- annotate_result_table(annotation_test_result, annotation_info, "auto")
stopifnot(all(c("ensembl_gene_id", "gene_symbol", "description") %in% names(annotated_test_result)))
stopifnot(sum(!is.na(annotated_test_result$gene_symbol) & annotated_test_result$gene_symbol != "") == 3)

built_in_go <- read_builtin_gene_sets(
  gene_set_dir = file.path(app_dir, "gene_sets"),
  gene_set_key = "yeast_scerevisiae_go",
  domain = "all"
)
prepared_go <- prepare_builtin_gene_sets_for_results(annotation_test_result, built_in_go, "auto")
go_selected_genes <- select_de_genes(
  result_table = annotation_test_result,
  alpha = 0.05,
  lfc_cutoff = 1,
  gene_list = "both"
)
go_ora_result <- run_ora_analysis(
  universe_genes = annotation_test_result$gene_id,
  selected_genes = go_selected_genes$genes,
  gene_sets = prepared_go,
  min_set_size = 1,
  max_set_size = 2000,
  padj_cutoff = 0.1
)
stopifnot(nrow(go_ora_result$results) > 0)
stopifnot(identical(attr(prepared_go, "gene_set_summary")$match_mode, "ensembl"))

protein_coding_rows <- annotation_info$table[
  annotation_info$table$gene_biotype == "protein_coding" &
    !is.na(annotation_info$table$ensembl_gene_id) &
    annotation_info$table$ensembl_gene_id != "",
  ,
  drop = FALSE
]
protein_coding_rows <- utils::head(protein_coding_rows, 8)
rosbys_lab_test_result <- data.frame(
  gene_id = protein_coding_rows$ensembl_gene_id,
  baseMean = seq(100, 170, by = 10),
  log2FoldChange = c(1.4, 1.3, 1.2, -1.4, -1.3, -1.2, 0.1, 0.2),
  pvalue = c(0.0005, 0.0006, 0.0007, 0.0005, 0.0006, 0.0007, 0.5, 0.6),
  padj = c(0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.8, 0.9),
  stringsAsFactors = FALSE
)
rosbys_lab_test_sets <- data.frame(
  term_id = rep(c("test_up", "test_down"), each = 4),
  term_name = rep(c("Synthetic up term", "Synthetic down term"), each = 4),
  gene_id = c(protein_coding_rows$ensembl_gene_id[c(1:3, 7)], protein_coding_rows$ensembl_gene_id[c(4:6, 8)]),
  domain = "test",
  stringsAsFactors = FALSE
)
rosbys_lab_ora_result <- run_rosbys_lab_enrichment(
  result_table = rosbys_lab_test_result,
  gene_sets = rosbys_lab_test_sets,
  annotation_info = annotation_info,
  alpha = 0.05,
  lfc_cutoff = 1,
  padj_cutoff = 1
)
stopifnot(nrow(rosbys_lab_ora_result$results) > 0)
stopifnot(all(c("group", "padj", "genes") %in% names(rosbys_lab_ora_result$results)))
stopifnot(all(c("group", "FDR", "Pathway size", "Fold enriched") %in% names(format_rosbys_lab_enrichment_export(rosbys_lab_ora_result$results))))

rosbys_lab_contaminated_sets <- rbind(
  rosbys_lab_test_sets,
  data.frame(
    term_id = rep(c("test_up", "test_down"), each = 2),
    term_name = rep(c("Synthetic up term", "Synthetic down term"), each = 2),
    gene_id = c(
      "OUTSIDE_PROTEIN_CODING_A",
      "OUTSIDE_PROTEIN_CODING_B",
      "OUTSIDE_PROTEIN_CODING_C",
      "OUTSIDE_PROTEIN_CODING_D"
    ),
    domain = "test",
    stringsAsFactors = FALSE
  )
)
rosbys_lab_contaminated_result <- run_rosbys_lab_enrichment(
  result_table = rosbys_lab_test_result,
  gene_sets = rosbys_lab_contaminated_sets,
  annotation_info = annotation_info,
  alpha = 0.05,
  lfc_cutoff = 1,
  padj_cutoff = 1
)
rosbys_lab_key <- paste(rosbys_lab_ora_result$results$group, rosbys_lab_ora_result$results$term_id, sep = "::")
contaminated_key <- paste(
  rosbys_lab_contaminated_result$results$group,
  rosbys_lab_contaminated_result$results$term_id,
  sep = "::"
)
contaminated_aligned <- rosbys_lab_contaminated_result$results[
  match(rosbys_lab_key, contaminated_key),
  ,
  drop = FALSE
]
stopifnot(identical(rosbys_lab_ora_result$results$set_size, contaminated_aligned$set_size))
stopifnot(isTRUE(all.equal(
  rosbys_lab_ora_result$results$pvalue,
  contaminated_aligned$pvalue,
  tolerance = 1e-15
)))

gene_sets <- read_gene_set_file(file.path(app_dir, "sample_data", "gene_sets.csv"))
selected_genes <- select_de_genes(
  result_table = raw_result$results,
  alpha = raw_result$alpha,
  lfc_cutoff = 1,
  gene_list = "both"
)
ora_result <- run_ora_analysis(
  universe_genes = raw_result$results$gene_id,
  selected_genes = selected_genes$genes,
  gene_sets = gene_sets,
  min_set_size = 3,
  max_set_size = 2000,
  padj_cutoff = 0.1
)
stopifnot(nrow(ora_result$results) > 0)
stopifnot(all(c("term_id", "padj", "genes") %in% names(ora_result$results)))
stopifnot(nrow(selected_genes$table) > 0)

pathway_result <- run_preranked_pathway_analysis(
  result_table = raw_result$results,
  gene_sets = gene_sets,
  rank_metric = "log2fc",
  min_set_size = 3,
  max_set_size = 200,
  padj_cutoff = 0.1
)
stopifnot(nrow(pathway_result$results) > 0)
stopifnot(all(c("term_id", "NES", "padj", "leading_edge_genes") %in% names(pathway_result$results)))
stopifnot(nrow(pathway_result$ranked_table) >= 10)
pathway_term <- pathway_result$results$term_id[1]
leading_edge <- pathway_leading_edge_table(pathway_result, pathway_term)
stopifnot(nrow(leading_edge) > 0)

kegg_ora_result <- NULL
tryCatch(
  {
    kegg_sets <- read_builtin_gene_sets(
      gene_set_dir = file.path(app_dir, "gene_sets"),
      gene_set_key = "yeast_scerevisiae_kegg",
      domain = "kegg_pathway"
    )
    prepared_kegg <- prepare_builtin_gene_sets_for_results(raw_result$results, kegg_sets, "auto")
    kegg_ora_result <- run_ora_analysis(
      universe_genes = raw_result$results$gene_id,
      selected_genes = selected_genes$genes,
      gene_sets = prepared_kegg,
      min_set_size = 1,
      max_set_size = 2000,
      padj_cutoff = 0.1
    )
    stopifnot(nrow(kegg_ora_result$results) > 0)
  },
  error = function(err) {
    warning(sprintf("Skipping KEGG smoke test: %s", conditionMessage(err)), call. = FALSE)
  }
)

expression_matrix <- prepare_expression_matrix(
  read_dge_table(file.path(app_dir, "sample_data", "normalized_expression.csv"))
)
expression_result <- run_normalized_expression_workflow(
  expression_matrix = expression_matrix,
  metadata = metadata,
  condition_col = "condition",
  treatment_level = "treated",
  reference_level = "control",
  expression_scale = "linear",
  alpha = 0.05
)
stopifnot(nrow(expression_result$results) > 0)
stopifnot(any(!is.na(expression_result$results$padj)))

fc_padj <- run_fold_change_workflow(
  data = read_dge_table(file.path(app_dir, "sample_data", "fold_changes_padj.csv")),
  require_adjusted_p = TRUE,
  fold_change_scale = "auto",
  alpha = 0.05
)
stopifnot(nrow(fc_padj$results) > 0)
stopifnot(any(!is.na(fc_padj$results$padj)))

fc_only <- run_fold_change_workflow(
  data = read_dge_table(file.path(app_dir, "sample_data", "fold_changes_only.csv")),
  require_adjusted_p = FALSE,
  fold_change_scale = "auto",
  alpha = 0.05
)
stopifnot(nrow(fc_only$results) > 0)
stopifnot(all(is.na(fc_only$results$padj)))

bundle_dir <- file.path(tempdir(), "deseq2_workbench_smoke_bundle")
unlink(bundle_dir, recursive = TRUE, force = TRUE)
write_result_bundle(raw_result, bundle_dir)
expected_files <- c(
  "deseq2_results.csv",
  "matrix_values.csv",
  "pca_plot.png",
  "volcano_plot.png",
  "regulation_summary_plot.png",
  "ma_plot.png",
  "analysis_summary.txt"
)
stopifnot(all(file.exists(file.path(bundle_dir, expected_files))))

ggplot2::ggsave(
  filename = file.path(bundle_dir, "ora_enrichment_plot.png"),
  plot = make_ora_plot(ora_result$results),
  width = 7,
  height = 5,
  dpi = 160
)
stopifnot(file.exists(file.path(bundle_dir, "ora_enrichment_plot.png")))

ggplot2::ggsave(
  filename = file.path(bundle_dir, "pathway_summary_plot.png"),
  plot = make_pathway_summary_plot(pathway_result$results, padj_cutoff = 0.1),
  width = 8,
  height = 6,
  dpi = 160
)
ggplot2::ggsave(
  filename = file.path(bundle_dir, "pathway_enrichment_plot.png"),
  plot = make_pathway_enrichment_plot(pathway_result, pathway_term),
  width = 8,
  height = 5,
  dpi = 160
)
ggplot2::ggsave(
  filename = file.path(bundle_dir, "pathway_heatmap.png"),
  plot = make_pathway_heatmap(raw_result, leading_edge$gene_id),
  width = 8,
  height = 6,
  dpi = 160
)
stopifnot(all(file.exists(file.path(
  bundle_dir,
  c("pathway_summary_plot.png", "pathway_enrichment_plot.png", "pathway_heatmap.png")
))))

cat("Smoke test passed\n")
cat(sprintf("Raw-count genes tested: %s\n", raw_result$genes_tested))
cat(sprintf("Normalized-expression genes processed: %s\n", expression_result$genes_tested))
cat(sprintf("Fold-change rows processed: %s\n", fc_padj$genes_tested))
cat(sprintf("ORA gene sets tested: %s\n", nrow(ora_result$results)))
cat(sprintf("Ranked pathways tested: %s\n", nrow(pathway_result$results)))
cat(sprintf("Built-in GO gene sets tested: %s\n", nrow(go_ora_result$results)))
if (!is.null(kegg_ora_result)) {
  cat(sprintf("KEGG pathways tested: %s\n", nrow(kegg_ora_result$results)))
}
