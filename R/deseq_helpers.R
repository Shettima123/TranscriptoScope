read_dge_table <- function(path, name = basename(path)) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Input file could not be found.", call. = FALSE)
  }

  lower_name <- tolower(name)
  sep <- if (grepl("\\.(tsv|txt)$", lower_name)) "\t" else ","

  data <- tryCatch(
    utils::read.table(
      path,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA", "NaN")
    ),
    error = function(err) {
      stop(sprintf("Could not read %s: %s", name, conditionMessage(err)), call. = FALSE)
    }
  )

  names(data) <- trimws(names(data))
  if (ncol(data) < 2) {
    stop(sprintf("%s must contain at least two columns.", name), call. = FALSE)
  }

  data
}

prepare_count_matrix <- function(counts_df) {
  if (ncol(counts_df) < 2) {
    stop("The count matrix needs a gene ID column plus one or more sample columns.", call. = FALSE)
  }

  gene_ids <- trimws(as.character(counts_df[[1]]))
  if (any(is.na(gene_ids) | gene_ids == "")) {
    stop("The first count-matrix column contains empty gene IDs.", call. = FALSE)
  }
  if (anyDuplicated(gene_ids)) {
    duplicated_gene <- gene_ids[duplicated(gene_ids)][1]
    stop(sprintf("Gene IDs must be unique. First duplicate: %s", duplicated_gene), call. = FALSE)
  }

  count_values <- counts_df[-1]
  numeric_values <- lapply(count_values, function(column) {
    if (is.factor(column)) {
      column <- as.character(column)
    }
    if (is.character(column)) {
      column <- gsub(",", "", column, fixed = TRUE)
    }
    suppressWarnings(as.numeric(column))
  })
  numeric_column_has_values <- vapply(numeric_values, function(column) any(!is.na(column)), logical(1))
  count_values <- count_values[numeric_column_has_values]
  numeric_values <- numeric_values[numeric_column_has_values]

  if (length(count_values) == 0) {
    stop("The count matrix does not contain numeric sample columns.", call. = FALSE)
  }

  sample_names <- trimws(names(count_values))
  if (any(is.na(sample_names) | sample_names == "")) {
    stop("All count-matrix sample columns need names.", call. = FALSE)
  }
  if (anyDuplicated(sample_names)) {
    duplicated_sample <- sample_names[duplicated(sample_names)][1]
    stop(sprintf("Sample names must be unique. First duplicate: %s", duplicated_sample), call. = FALSE)
  }
  names(count_values) <- sample_names

  count_matrix <- as.matrix(as.data.frame(numeric_values, check.names = FALSE))
  rownames(count_matrix) <- gene_ids
  colnames(count_matrix) <- sample_names

  if (anyNA(count_matrix)) {
    stop("The count matrix contains missing or non-numeric count values.", call. = FALSE)
  }
  if (any(count_matrix < 0)) {
    stop("The count matrix contains negative values. DESeq2 expects non-negative raw counts.", call. = FALSE)
  }
  if (any(abs(count_matrix - round(count_matrix)) > 1e-6)) {
    stop(
      "The count matrix contains decimal values. DESeq2 expects raw integer counts, not TPM/FPKM/RPKM.",
      call. = FALSE
    )
  }
  if (!any(rowSums(count_matrix) > 0)) {
    stop("All genes have zero counts.", call. = FALSE)
  }

  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  count_matrix
}

prepare_expression_matrix <- function(expression_df) {
  if (ncol(expression_df) < 2) {
    stop("The expression matrix needs a gene ID column plus one or more sample columns.", call. = FALSE)
  }

  gene_ids <- trimws(as.character(expression_df[[1]]))
  if (any(is.na(gene_ids) | gene_ids == "")) {
    stop("The first expression-matrix column contains empty gene IDs.", call. = FALSE)
  }
  if (anyDuplicated(gene_ids)) {
    duplicated_gene <- gene_ids[duplicated(gene_ids)][1]
    stop(sprintf("Gene IDs must be unique. First duplicate: %s", duplicated_gene), call. = FALSE)
  }

  expression_values <- expression_df[-1]
  numeric_values <- lapply(expression_values, function(column) {
    if (is.factor(column)) {
      column <- as.character(column)
    }
    if (is.character(column)) {
      column <- gsub(",", "", column, fixed = TRUE)
    }
    suppressWarnings(as.numeric(column))
  })
  numeric_column_has_values <- vapply(numeric_values, function(column) any(!is.na(column)), logical(1))
  expression_values <- expression_values[numeric_column_has_values]
  numeric_values <- numeric_values[numeric_column_has_values]

  if (length(expression_values) == 0) {
    stop("The expression matrix does not contain numeric sample columns.", call. = FALSE)
  }

  sample_names <- trimws(names(expression_values))
  if (any(is.na(sample_names) | sample_names == "")) {
    stop("All expression-matrix sample columns need names.", call. = FALSE)
  }
  if (anyDuplicated(sample_names)) {
    duplicated_sample <- sample_names[duplicated(sample_names)][1]
    stop(sprintf("Sample names must be unique. First duplicate: %s", duplicated_sample), call. = FALSE)
  }
  names(expression_values) <- sample_names

  expression_matrix <- as.matrix(as.data.frame(numeric_values, check.names = FALSE))
  rownames(expression_matrix) <- gene_ids
  colnames(expression_matrix) <- sample_names

  if (anyNA(expression_matrix)) {
    stop("The expression matrix contains missing or non-numeric values.", call. = FALSE)
  }
  if (nrow(expression_matrix) < 2) {
    stop("The expression matrix needs at least two genes.", call. = FALSE)
  }
  if (ncol(expression_matrix) < 2) {
    stop("The expression matrix needs at least two samples.", call. = FALSE)
  }

  expression_matrix
}

prepare_sample_metadata <- function(metadata_df) {
  if (ncol(metadata_df) < 2) {
    stop("The metadata table needs a sample ID column plus at least one analysis column.", call. = FALSE)
  }

  clean_names <- tolower(gsub("[^a-z0-9]", "", names(metadata_df)))
  sample_col <- match(TRUE, clean_names %in% c("sampleid", "sample", "samplename"))
  if (is.na(sample_col)) {
    sample_col <- 1
  }

  names(metadata_df) <- make.unique(trimws(names(metadata_df)), sep = "_")
  names(metadata_df)[sample_col] <- "sample_id"

  sample_ids <- trimws(as.character(metadata_df[[sample_col]]))
  if (any(is.na(sample_ids) | sample_ids == "")) {
    stop("The metadata sample ID column contains empty sample names.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    duplicated_sample <- sample_ids[duplicated(sample_ids)][1]
    stop(sprintf("Metadata sample IDs must be unique. First duplicate: %s", duplicated_sample), call. = FALSE)
  }

  metadata_df$sample_id <- sample_ids
  rownames(metadata_df) <- sample_ids
  metadata_df
}

infer_condition_from_sample_names <- function(sample_names) {
  labels <- trimws(sample_names)
  labels <- sub("([._ -]?(rep|replicate)[._ -]*[0-9]+)$", "", labels, ignore.case = TRUE)
  labels <- sub("([._ -]?r[0-9]+)$", "", labels, ignore.case = TRUE)
  labels <- sub("([._ -]?[0-9]+)$", "", labels, ignore.case = TRUE)
  labels <- gsub("[._ -]+$", "", labels)

  if (length(unique(labels)) < 2 && any(grepl("_", sample_names, fixed = TRUE))) {
    labels <- sub("_[^_]+$", "", sample_names)
  }
  labels[is.na(labels) | labels == ""] <- sample_names[is.na(labels) | labels == ""]
  labels
}

infer_sample_metadata <- function(sample_names) {
  data.frame(
    sample_id = sample_names,
    condition = infer_condition_from_sample_names(sample_names),
    stringsAsFactors = FALSE,
    row.names = sample_names
  )
}

align_metadata_to_counts <- function(count_matrix, metadata) {
  missing_metadata <- setdiff(colnames(count_matrix), rownames(metadata))
  extra_metadata <- setdiff(rownames(metadata), colnames(count_matrix))

  if (length(missing_metadata) > 0) {
    stop(
      sprintf(
        "Metadata is missing sample(s) found in the count matrix: %s",
        paste(missing_metadata, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  warnings <- character()
  if (length(extra_metadata) > 0) {
    warnings <- c(
      warnings,
      sprintf(
        "Metadata has extra sample(s) not present in the count matrix and they will be ignored: %s",
        paste(extra_metadata, collapse = ", ")
      )
    )
  }

  aligned <- metadata[colnames(count_matrix), , drop = FALSE]
  attr(aligned, "dge_warnings") <- warnings
  aligned
}

metadata_for_matrix <- function(matrix_values, metadata = NULL) {
  if (is.null(metadata)) {
    inferred <- infer_sample_metadata(colnames(matrix_values))
    attr(inferred, "dge_warnings") <- paste(
      "No metadata file was supplied, so sample groups were inferred from column names:",
      paste(unique(inferred$condition), collapse = ", ")
    )
    return(inferred)
  }

  tryCatch(
    align_metadata_to_counts(matrix_values, metadata),
    error = function(err) {
      inferred <- infer_sample_metadata(colnames(matrix_values))
      attr(inferred, "dge_warnings") <- paste(
        "The supplied metadata did not match the data-file sample columns, so it was ignored and groups were inferred from column names:",
        paste(unique(inferred$condition), collapse = ", ")
      )
      inferred
    }
  )
}

cpm_matrix <- function(count_matrix) {
  library_sizes <- colSums(count_matrix)
  if (any(library_sizes <= 0)) {
    stop("At least one sample has zero total counts, so CPM cannot be calculated.", call. = FALSE)
  }

  sweep(count_matrix, 2, library_sizes, "/") * 1e6
}

matrix_to_export_table <- function(matrix_values) {
  export_table <- as.data.frame(matrix_values, check.names = FALSE)
  export_table$gene_id <- rownames(export_table)
  export_table <- export_table[c("gene_id", setdiff(names(export_table), "gene_id"))]
  rownames(export_table) <- NULL
  export_table
}

read_annotation_manifest <- function(annotation_dir) {
  manifest_file <- file.path(annotation_dir, "manifest.csv")
  if (!file.exists(manifest_file)) {
    return(data.frame())
  }

  manifest <- utils::read.csv(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("key", "label", "scientific_name", "dataset", "taxonomy_id", "assembly", "file", "source", "rows")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0) {
    stop(sprintf("Annotation manifest is missing column(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  manifest
}

read_gene_annotation <- function(annotation_dir, annotation_key) {
  manifest <- read_annotation_manifest(annotation_dir)
  if (nrow(manifest) == 0) {
    stop("No built-in annotations are available.", call. = FALSE)
  }

  row <- manifest[manifest$key == annotation_key, , drop = FALSE]
  if (nrow(row) == 0) {
    stop("Choose a valid built-in annotation.", call. = FALSE)
  }

  annotation_file <- file.path(annotation_dir, row$file[1])
  if (!file.exists(annotation_file)) {
    stop(sprintf("Annotation file is missing: %s", row$file[1]), call. = FALSE)
  }

  annotation <- utils::read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("ensembl_gene_id", "gene_symbol", "description", "gene_biotype", "chromosome", "start", "end")
  missing <- setdiff(required, names(annotation))
  if (length(missing) > 0) {
    stop(sprintf("Annotation file is missing column(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  list(
    table = annotation,
    manifest = row
  )
}

read_gene_set_manifest <- function(gene_set_dir) {
  manifest_file <- file.path(gene_set_dir, "manifest.csv")
  if (!file.exists(manifest_file)) {
    return(data.frame())
  }

  manifest <- utils::read.csv(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("key", "annotation_key", "label", "dataset", "taxonomy_id", "assembly", "file", "rows", "terms")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0) {
    stop(sprintf("Gene set manifest is missing column(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }
  if (!"collection" %in% names(manifest)) {
    manifest$collection <- "go"
  }
  if (!"kegg_organism" %in% names(manifest)) {
    manifest$kegg_organism <- ""
  }
  manifest$collection[is.na(manifest$collection) | manifest$collection == ""] <- "go"
  manifest$kegg_organism[is.na(manifest$kegg_organism)] <- ""
  manifest
}

.dge_gene_set_cache <- new.env(parent = emptyenv())

download_kegg_lines <- function(url) {
  output_file <- tempfile(fileext = ".txt")
  on.exit(unlink(output_file), add = TRUE)

  status <- tryCatch(
    utils::download.file(url, output_file, quiet = TRUE, mode = "wb"),
    error = function(err) err
  )
  if (inherits(status, "error") || !identical(status, 0L)) {
    stop(
      sprintf("Could not download KEGG data from %s. Check internet access or try again later.", url),
      call. = FALSE
    )
  }

  lines <- readLines(output_file, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    stop(sprintf("KEGG returned no data from %s.", url), call. = FALSE)
  }
  lines
}

read_kegg_tsv <- function(url) {
  lines <- download_kegg_lines(url)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  max_cols <- max(vapply(split_lines, length, integer(1)))
  rows <- lapply(split_lines, function(parts) {
    length(parts) <- max_cols
    parts[is.na(parts)] <- ""
    parts
  })
  data.frame(do.call(rbind, rows), stringsAsFactors = FALSE, check.names = FALSE)
}

strip_kegg_prefix <- function(ids) {
  sub("^[^:]+:", "", ids)
}

extract_kegg_symbol <- function(descriptions) {
  symbols <- vapply(strsplit(descriptions, ";", fixed = TRUE), function(parts) {
    alias <- trimws(parts[1])
    alias <- trimws(strsplit(alias, ",", fixed = TRUE)[[1]][1])
    if (!nzchar(alias) || grepl("\\s", alias) || grepl("uncharacterized|hypothetical|putative protein", alias, ignore.case = TRUE)) {
      return("")
    }
    alias
  }, character(1))
  symbols
}

read_kegg_gene_sets <- function(gene_set_dir, manifest_row) {
  organism <- manifest_row$kegg_organism[1]
  if (is.na(organism) || !nzchar(organism)) {
    stop("The selected KEGG database is missing a KEGG organism code.", call. = FALSE)
  }

  cache_dir <- file.path(gene_set_dir, "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, sprintf("%s.csv", manifest_row$key[1]))
  if (file.exists(cache_file)) {
    gene_sets <- utils::read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE)
    return(gene_sets)
  }

  base_url <- "https://rest.kegg.jp"
  pathway_list <- read_kegg_tsv(sprintf("%s/list/pathway/%s", base_url, organism))
  pathway_links <- read_kegg_tsv(sprintf("%s/link/%s/pathway", base_url, organism))
  gene_list <- read_kegg_tsv(sprintf("%s/list/%s", base_url, organism))

  names(pathway_list)[1:2] <- c("term_id", "term_name")
  pathway_names <- pathway_list[, c("term_id", "term_name"), drop = FALSE]
  pathway_names$term_id <- sub("^path:", "", pathway_names$term_id)

  names(pathway_links)[1:2] <- c("term_id", "kegg_gene_id")
  pathway_links <- pathway_links[, c("term_id", "kegg_gene_id"), drop = FALSE]
  pathway_links$term_id <- sub("^path:", "", pathway_links$term_id)

  names(gene_list)[1] <- "kegg_gene_id"
  description_col <- if (ncol(gene_list) >= 4) 4 else if (ncol(gene_list) >= 2) 2 else 1
  gene_lookup <- data.frame(
    kegg_gene_id = gene_list$kegg_gene_id,
    gene_id = strip_kegg_prefix(gene_list$kegg_gene_id),
    gene_symbol = extract_kegg_symbol(gene_list[[description_col]]),
    stringsAsFactors = FALSE
  )
  gene_lookup <- gene_lookup[!duplicated(gene_lookup$kegg_gene_id), , drop = FALSE]

  gene_sets <- merge(pathway_links, pathway_names, by = "term_id", all.x = TRUE, sort = FALSE)
  gene_sets <- merge(gene_sets, gene_lookup, by = "kegg_gene_id", all.x = TRUE, sort = FALSE)
  gene_sets$term_name[is.na(gene_sets$term_name) | gene_sets$term_name == ""] <- gene_sets$term_id[is.na(gene_sets$term_name) | gene_sets$term_name == ""]
  gene_sets$gene_id[is.na(gene_sets$gene_id)] <- strip_kegg_prefix(gene_sets$kegg_gene_id[is.na(gene_sets$gene_id)])
  gene_sets$gene_symbol[is.na(gene_sets$gene_symbol)] <- ""
  gene_sets$domain <- "kegg_pathway"
  gene_sets <- unique(gene_sets[c("gene_id", "gene_symbol", "term_id", "term_name", "domain")])
  gene_sets <- gene_sets[!is.na(gene_sets$gene_id) & gene_sets$gene_id != "", , drop = FALSE]
  if (nrow(gene_sets) == 0) {
    stop("KEGG returned no pathway-to-gene mappings for the selected organism.", call. = FALSE)
  }

  utils::write.csv(gene_sets, cache_file, row.names = FALSE, na = "")
  gene_sets
}

read_builtin_gene_sets <- function(gene_set_dir, gene_set_key, domain = "all") {
  manifest <- read_gene_set_manifest(gene_set_dir)
  if (nrow(manifest) == 0) {
    stop("No built-in gene set databases are available.", call. = FALSE)
  }

  row <- manifest[manifest$key == gene_set_key, , drop = FALSE]
  if (nrow(row) == 0) {
    stop("Choose a valid built-in gene set database.", call. = FALSE)
  }

  cache_key <- paste(gene_set_key, domain, sep = "::")
  if (exists(cache_key, envir = .dge_gene_set_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .dge_gene_set_cache, inherits = FALSE))
  }

  if (identical(row$collection[1], "kegg")) {
    gene_sets <- read_kegg_gene_sets(gene_set_dir, row)
  } else {
    gene_set_file <- file.path(gene_set_dir, row$file[1])
    if (!file.exists(gene_set_file)) {
      stop(sprintf("Built-in gene set file is missing: %s", row$file[1]), call. = FALSE)
    }

    gene_sets <- utils::read.csv(
      gene_set_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  required <- c("gene_id", "gene_symbol", "term_id", "term_name", "domain")
  missing <- setdiff(required, names(gene_sets))
  if (length(missing) > 0) {
    stop(sprintf("Built-in gene set file is missing column(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  if (!identical(domain, "all") && !identical(row$collection[1], "kegg")) {
    gene_sets <- gene_sets[gene_sets$domain == domain, , drop = FALSE]
  }
  gene_sets <- gene_sets[!is.na(gene_sets$term_id) & gene_sets$term_id != "", , drop = FALSE]
  if (nrow(gene_sets) == 0) {
    stop("No built-in gene sets remain after domain filtering.", call. = FALSE)
  }

  result <- list(
    table = gene_sets,
    manifest = row,
    domain = domain
  )
  assign(cache_key, result, envir = .dge_gene_set_cache)
  result
}

annotation_match_stats <- function(result_table, annotation) {
  result_ids <- normalize_gene_ids(result_table$gene_id, case_sensitive = FALSE)
  ensembl_ids <- normalize_gene_ids(annotation$ensembl_gene_id, case_sensitive = FALSE)
  symbols <- normalize_gene_ids(annotation$gene_symbol, case_sensitive = FALSE)
  symbols <- symbols[nzchar(symbols)]

  c(
    ensembl = sum(result_ids %in% ensembl_ids),
    symbol = sum(result_ids %in% symbols)
  )
}

gene_set_match_stats <- function(result_table, gene_sets) {
  result_ids <- normalize_gene_ids(result_table$gene_id, case_sensitive = FALSE)
  gene_ids <- normalize_gene_ids(gene_sets$gene_id, case_sensitive = FALSE)
  symbols <- normalize_gene_ids(gene_sets$gene_symbol, case_sensitive = FALSE)
  symbols <- symbols[nzchar(symbols)]

  c(
    ensembl = sum(result_ids %in% gene_ids),
    symbol = sum(result_ids %in% symbols)
  )
}

prepare_builtin_gene_sets_for_results <- function(result_table, built_in_gene_sets, match_mode = "auto") {
  gene_sets <- built_in_gene_sets$table
  resolved_mode <- match_mode
  if (identical(resolved_mode, "auto")) {
    stats <- gene_set_match_stats(result_table, gene_sets)
    resolved_mode <- if (stats[["symbol"]] > stats[["ensembl"]]) "symbol" else "ensembl"
  }

  gene_col <- if (identical(resolved_mode, "symbol")) "gene_symbol" else "gene_id"
  prepared <- data.frame(
    term_id = gene_sets$term_id,
    term_name = ifelse(
      !is.na(gene_sets$term_name) & gene_sets$term_name != "",
      gene_sets$term_name,
      gene_sets$term_id
    ),
    gene_id = gene_sets[[gene_col]],
    domain = gene_sets$domain,
    stringsAsFactors = FALSE
  )
  prepared$gene_id <- trimws(as.character(prepared$gene_id))
  prepared <- prepared[!is.na(prepared$gene_id) & prepared$gene_id != "", , drop = FALSE]
  prepared <- unique(prepared)

  stats <- gene_set_match_stats(result_table, gene_sets)
  attr(prepared, "gene_set_summary") <- list(
    label = built_in_gene_sets$manifest$label[1],
    dataset = built_in_gene_sets$manifest$dataset[1],
    assembly = built_in_gene_sets$manifest$assembly[1],
    domain = built_in_gene_sets$domain,
    match_mode = resolved_mode,
    matched = stats[[resolved_mode]],
    total = nrow(result_table),
    terms = length(unique(prepared$term_id)),
    rows = nrow(prepared)
  )
  prepared
}

resolve_annotation_match_mode <- function(result_table, annotation, match_mode = "auto") {
  if (identical(match_mode, "ensembl") || identical(match_mode, "symbol")) {
    return(match_mode)
  }
  stats <- annotation_match_stats(result_table, annotation)
  if (stats[["symbol"]] > stats[["ensembl"]]) "symbol" else "ensembl"
}

annotate_result_table <- function(result_table, annotation_info = NULL, match_mode = "auto") {
  if (is.null(annotation_info)) {
    return(result_table)
  }

  annotation <- annotation_info$table
  resolved_mode <- resolve_annotation_match_mode(result_table, annotation, match_mode)
  annotation_key_col <- if (identical(resolved_mode, "symbol")) "gene_symbol" else "ensembl_gene_id"

  annotation_key <- normalize_gene_ids(annotation[[annotation_key_col]], case_sensitive = FALSE)
  keep <- !is.na(annotation_key) & annotation_key != "" & !duplicated(annotation_key)
  annotation_keyed <- annotation[keep, , drop = FALSE]
  annotation_keyed$.annotation_key <- annotation_key[keep]

  result <- result_table
  result$.annotation_key <- normalize_gene_ids(result$gene_id, case_sensitive = FALSE)
  result$.row_order <- seq_len(nrow(result))

  annotated <- merge(
    result,
    annotation_keyed,
    by = ".annotation_key",
    all.x = TRUE,
    sort = FALSE
  )
  annotated <- annotated[order(annotated$.row_order), , drop = FALSE]

  annotation_cols <- c("ensembl_gene_id", "gene_symbol", "description", "gene_biotype", "chromosome", "start", "end")
  result_cols <- setdiff(names(result_table), annotation_cols)
  annotated <- annotated[c(result_cols, annotation_cols)]
  rownames(annotated) <- NULL

  matched <- sum(!is.na(annotated[[annotation_key_col]]) & annotated[[annotation_key_col]] != "")
  attr(annotated, "annotation_summary") <- list(
    label = annotation_info$manifest$label[1],
    dataset = annotation_info$manifest$dataset[1],
    assembly = annotation_info$manifest$assembly[1],
    match_mode = resolved_mode,
    matched = matched,
    total = nrow(result_table)
  )
  annotated
}

preprocess_count_matrix <- function(
  count_matrix,
  metadata = NULL,
  min_cpm = 0.5,
  min_samples = 1,
  pseudocount = 4
) {
  if (min_cpm < 0) {
    stop("Minimum CPM must be 0 or greater.", call. = FALSE)
  }
  if (min_samples < 1) {
    stop("Minimum samples must be at least 1.", call. = FALSE)
  }
  if (pseudocount < 0) {
    stop("Pseudocount must be 0 or greater.", call. = FALSE)
  }

  min_samples <- min(as.integer(round(min_samples)), ncol(count_matrix))
  aligned_metadata <- metadata_for_matrix(count_matrix, metadata)
  raw_totals <- colSums(count_matrix)
  cpm <- cpm_matrix(count_matrix)
  keep <- rowSums(cpm >= min_cpm) >= min_samples

  if (!any(keep)) {
    stop("No genes pass the preprocessing filter. Lower the CPM threshold or minimum sample count.", call. = FALSE)
  }

  filtered_counts <- count_matrix[keep, , drop = FALSE]
  transformed <- log2(cpm_matrix(filtered_counts) + pseudocount)
  filtered_totals <- colSums(filtered_counts)

  sample_summary <- data.frame(
    sample_id = colnames(count_matrix),
    total_counts = raw_totals,
    total_counts_after_filter = filtered_totals,
    counts_removed = raw_totals - filtered_totals,
    stringsAsFactors = FALSE
  )

  list(
    filtered_counts = filtered_counts,
    transformed_expression = transformed,
    metadata = aligned_metadata,
    sample_summary = sample_summary,
    genes_before = nrow(count_matrix),
    genes_after = nrow(filtered_counts),
    genes_removed = nrow(count_matrix) - nrow(filtered_counts),
    min_cpm = min_cpm,
    min_samples = min_samples,
    pseudocount = pseudocount,
    warnings = attr(aligned_metadata, "dge_warnings")
  )
}

analysis_columns <- function(metadata) {
  setdiff(names(metadata), "sample_id")
}

condition_levels <- function(metadata, condition_col) {
  values <- trimws(as.character(metadata[[condition_col]]))
  values <- values[!is.na(values) & values != ""]
  unique(values)
}

validate_analysis_settings <- function(
  count_matrix,
  metadata,
  condition_col,
  treatment_level,
  reference_level,
  batch_col = NULL,
  min_total_count = 10
) {
  warnings <- character()

  if (!condition_col %in% names(metadata)) {
    stop("Choose a valid condition column.", call. = FALSE)
  }

  condition <- trimws(as.character(metadata[[condition_col]]))
  if (any(is.na(condition) | condition == "")) {
    stop("The selected condition column contains missing values.", call. = FALSE)
  }

  levels_found <- unique(condition)
  if (length(levels_found) < 2) {
    stop("The selected condition column needs at least two groups.", call. = FALSE)
  }
  if (!reference_level %in% levels_found) {
    stop("The reference group was not found in the selected condition column.", call. = FALSE)
  }
  if (!treatment_level %in% levels_found) {
    stop("The comparison group was not found in the selected condition column.", call. = FALSE)
  }
  if (identical(reference_level, treatment_level)) {
    stop("Reference and comparison groups must be different.", call. = FALSE)
  }

  replicate_counts <- table(condition)
  if (any(replicate_counts < 2)) {
    warnings <- c(
      warnings,
      "At least one condition has fewer than two replicates. DESeq2 can run, but interpretation will be fragile."
    )
  }

  if (!is.null(batch_col) && nzchar(batch_col)) {
    if (!batch_col %in% names(metadata)) {
      stop("Choose a valid batch column.", call. = FALSE)
    }
    batch <- trimws(as.character(metadata[[batch_col]]))
    if (any(is.na(batch) | batch == "")) {
      stop("The selected batch column contains missing values.", call. = FALSE)
    }
    if (length(unique(batch)) < 2) {
      warnings <- c(warnings, "The selected batch column has only one level, so it does not add information.")
    }
  }

  genes_kept <- sum(rowSums(count_matrix) >= min_total_count)
  if (genes_kept < 2) {
    stop(
      sprintf("Fewer than two genes pass the minimum total count filter of %s.", min_total_count),
      call. = FALSE
    )
  }

  warnings
}

run_deseq2_workflow <- function(
  count_matrix,
  metadata,
  condition_col,
  treatment_level,
  reference_level,
  batch_col = NULL,
  min_total_count = 10,
  alpha = 0.05
) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("S4Vectors is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }

  batch_col <- if (is.null(batch_col) || identical(batch_col, "__none__")) NULL else batch_col
  input_warnings <- attr(metadata, "dge_warnings")
  aligned_metadata <- align_metadata_to_counts(count_matrix, metadata)
  workflow_warnings <- c(input_warnings, attr(aligned_metadata, "dge_warnings"))
  workflow_warnings <- c(
    workflow_warnings,
    validate_analysis_settings(
      count_matrix = count_matrix,
      metadata = aligned_metadata,
      condition_col = condition_col,
      treatment_level = treatment_level,
      reference_level = reference_level,
      batch_col = batch_col,
      min_total_count = min_total_count
    )
  )

  keep <- rowSums(count_matrix) >= min_total_count
  filtered_counts <- count_matrix[keep, , drop = FALSE]

  condition <- factor(trimws(as.character(aligned_metadata[[condition_col]])))
  condition <- stats::relevel(condition, ref = reference_level)
  analysis_metadata <- data.frame(condition = condition, row.names = rownames(aligned_metadata))

  if (!is.null(batch_col)) {
    analysis_metadata$batch <- factor(trimws(as.character(aligned_metadata[[batch_col]])))
    design <- stats::as.formula("~ batch + condition")
  } else {
    design <- stats::as.formula("~ condition")
  }

  model_matrix <- stats::model.matrix(design, data = analysis_metadata)
  if (qr(model_matrix)$rank < ncol(model_matrix)) {
    stop(
      "The selected design cannot be fit because condition and batch are confounded. Choose a simpler design or different batch column.",
      call. = FALSE
    )
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = analysis_metadata,
    design = design
  )
  dds <- tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(err) {
      message <- conditionMessage(err)
      if (!grepl("dispersion estimates", message)) {
        stop(err)
      }

      workflow_warnings <<- c(
        workflow_warnings,
        "DESeq2's default dispersion fit failed, so the workflow used gene-wise dispersion estimates."
      )
      dds <- DESeq2::estimateSizeFactors(dds)
      dds <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
      dispersion_setter <- get("dispersions<-", envir = asNamespace("DESeq2"))
      dds <- dispersion_setter(dds, value = S4Vectors::mcols(dds)$dispGeneEst)
      DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    }
  )

  deseq_result <- DESeq2::results(
    dds,
    contrast = c("condition", treatment_level, reference_level),
    alpha = alpha
  )

  result_table <- as.data.frame(deseq_result)
  result_table$gene_id <- rownames(result_table)
  result_table <- result_table[c("gene_id", setdiff(names(result_table), "gene_id"))]
  result_table <- result_table[order(result_table$padj, result_table$pvalue, na.last = TRUE), , drop = FALSE]
  rownames(result_table) <- NULL

  normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE), check.names = FALSE)
  normalized_counts$gene_id <- rownames(normalized_counts)
  normalized_counts <- normalized_counts[c("gene_id", setdiff(names(normalized_counts), "gene_id"))]
  rownames(normalized_counts) <- NULL

  vst_data <- tryCatch(
    {
      if (nrow(dds) < 1000) {
        workflow_warnings <- c(
          workflow_warnings,
          "The dataset has fewer than 1,000 genes after filtering, so PCA used varianceStabilizingTransformation instead of vst."
        )
        DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
      } else {
        DESeq2::vst(dds, blind = FALSE)
      }
    },
    error = function(err) {
      workflow_warnings <<- c(
        workflow_warnings,
        sprintf("PCA used log2 normalized counts because the DESeq2 transform failed: %s", conditionMessage(err))
      )
      NULL
    }
  )

  if (is.null(vst_data)) {
    pca_matrix <- log2(as.matrix(normalized_counts[setdiff(names(normalized_counts), "gene_id")]) + 1)
    rownames(pca_matrix) <- normalized_counts$gene_id
    pca <- stats::prcomp(t(pca_matrix), center = TRUE, scale. = FALSE)
    explained <- pca$sdev^2 / sum(pca$sdev^2)
    percent_var <- round(100 * explained[1:2])
    pca_data <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      condition = analysis_metadata$condition,
      name = rownames(analysis_metadata),
      stringsAsFactors = FALSE
    )
  } else {
    pca_data <- DESeq2::plotPCA(vst_data, intgroup = "condition", returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))
  }

  list(
    dds = dds,
    deseq_result = deseq_result,
    results = result_table,
    normalized_counts = normalized_counts,
    pca_data = pca_data,
    percent_var = percent_var,
    warnings = unique(workflow_warnings),
    genes_tested = nrow(filtered_counts),
    samples_tested = ncol(filtered_counts),
    condition_col = condition_col,
    treatment_level = treatment_level,
    reference_level = reference_level,
    batch_col = batch_col,
    alpha = alpha,
    input_type = "counts",
    workflow_label = "DESeq2"
  )
}

resolve_expression_scale <- function(expression_matrix, expression_scale = "auto") {
  if (identical(expression_scale, "log2") || identical(expression_scale, "linear")) {
    return(expression_scale)
  }

  finite_values <- expression_matrix[is.finite(expression_matrix)]
  if (length(finite_values) == 0) {
    return("linear")
  }
  if (any(finite_values < 0)) {
    return("log2")
  }
  if (stats::median(finite_values) <= 20 && stats::quantile(finite_values, 0.99, names = FALSE) <= 100) {
    return("log2")
  }
  "linear"
}

make_pca_data_from_matrix <- function(matrix_values, metadata, condition_col) {
  if (nrow(matrix_values) < 2 || ncol(matrix_values) < 2) {
    return(NULL)
  }

  pca <- stats::prcomp(t(matrix_values), center = TRUE, scale. = FALSE)
  explained <- pca$sdev^2 / sum(pca$sdev^2)
  pca_scores <- pca$x

  if (ncol(pca_scores) < 2) {
    pca_scores <- cbind(pca_scores[, 1], 0)
  }

  pca_data <- data.frame(
    PC1 = pca_scores[, 1],
    PC2 = pca_scores[, 2],
    condition = metadata[[condition_col]],
    name = rownames(metadata),
    stringsAsFactors = FALSE
  )

  list(
    pca_data = pca_data,
    percent_var = round(100 * explained[1:2])
  )
}

safe_row_t_test <- function(treatment_values, reference_values) {
  if (length(treatment_values) < 2 || length(reference_values) < 2) {
    return(NA_real_)
  }

  if (stats::sd(treatment_values) == 0 && stats::sd(reference_values) == 0) {
    return(if (mean(treatment_values) == mean(reference_values)) 1 else NA_real_)
  }

  tryCatch(
    stats::t.test(treatment_values, reference_values)$p.value,
    error = function(err) NA_real_
  )
}

sort_result_table <- function(result_table) {
  if ("padj" %in% names(result_table) && any(!is.na(result_table$padj))) {
    result_table <- result_table[order(result_table$padj, result_table$pvalue, na.last = TRUE), , drop = FALSE]
  } else if ("log2FoldChange" %in% names(result_table)) {
    result_table <- result_table[order(abs(result_table$log2FoldChange), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  rownames(result_table) <- NULL
  result_table
}

run_normalized_expression_workflow <- function(
  expression_matrix,
  metadata,
  condition_col,
  treatment_level,
  reference_level,
  expression_scale = "auto",
  alpha = 0.05,
  pseudocount = 1
) {
  aligned_metadata <- metadata_for_matrix(expression_matrix, metadata)
  workflow_warnings <- c(attr(metadata, "dge_warnings"), attr(aligned_metadata, "dge_warnings"))
  workflow_warnings <- c(
    workflow_warnings,
    validate_analysis_settings(
      count_matrix = abs(round(expression_matrix)),
      metadata = aligned_metadata,
      condition_col = condition_col,
      treatment_level = treatment_level,
      reference_level = reference_level,
      batch_col = NULL,
      min_total_count = 0
    ),
    "Normalized expression mode does not run DESeq2. It reports group means, log2 fold change, and Welch t-test p-values when replicates are available."
  )

  resolved_scale <- resolve_expression_scale(expression_matrix, expression_scale)
  condition <- trimws(as.character(aligned_metadata[[condition_col]]))
  reference_samples <- rownames(aligned_metadata)[condition == reference_level]
  treatment_samples <- rownames(aligned_metadata)[condition == treatment_level]

  reference_values <- expression_matrix[, reference_samples, drop = FALSE]
  treatment_values <- expression_matrix[, treatment_samples, drop = FALSE]
  reference_mean <- rowMeans(reference_values)
  treatment_mean <- rowMeans(treatment_values)

  if (identical(resolved_scale, "linear")) {
    if (any(expression_matrix < 0)) {
      stop("Linear normalized expression values cannot be negative. Choose log2-transformed expression scale.", call. = FALSE)
    }
    log2_fold_change <- log2((treatment_mean + pseudocount) / (reference_mean + pseudocount))
    pca_matrix <- log2(expression_matrix + pseudocount)
  } else {
    log2_fold_change <- treatment_mean - reference_mean
    pca_matrix <- expression_matrix
  }

  pvalues <- vapply(
    seq_len(nrow(expression_matrix)),
    function(row_index) safe_row_t_test(treatment_values[row_index, ], reference_values[row_index, ]),
    FUN.VALUE = numeric(1)
  )
  padj <- stats::p.adjust(pvalues, method = "BH")

  result_table <- data.frame(
    gene_id = rownames(expression_matrix),
    baseMean = rowMeans(expression_matrix),
    log2FoldChange = log2_fold_change,
    lfcSE = NA_real_,
    stat = NA_real_,
    pvalue = pvalues,
    padj = padj,
    referenceMean = reference_mean,
    treatmentMean = treatment_mean,
    stringsAsFactors = FALSE
  )
  result_table <- sort_result_table(result_table)

  matrix_export <- as.data.frame(expression_matrix, check.names = FALSE)
  matrix_export$gene_id <- rownames(matrix_export)
  matrix_export <- matrix_export[c("gene_id", setdiff(names(matrix_export), "gene_id"))]
  rownames(matrix_export) <- NULL

  pca <- make_pca_data_from_matrix(pca_matrix, aligned_metadata, condition_col)

  list(
    dds = NULL,
    deseq_result = NULL,
    results = result_table,
    normalized_counts = matrix_export,
    pca_data = if (is.null(pca)) NULL else pca$pca_data,
    percent_var = if (is.null(pca)) NULL else pca$percent_var,
    warnings = unique(c(workflow_warnings, sprintf("Expression scale used: %s.", resolved_scale))),
    genes_tested = nrow(expression_matrix),
    samples_tested = ncol(expression_matrix),
    condition_col = condition_col,
    treatment_level = treatment_level,
    reference_level = reference_level,
    batch_col = NULL,
    alpha = alpha,
    input_type = "normalized",
    workflow_label = "Normalized expression"
  )
}

clean_column_name <- function(name) {
  tolower(gsub("[^a-z0-9]", "", name))
}

find_named_column <- function(data, candidates) {
  cleaned <- clean_column_name(names(data))
  matched <- match(candidates, cleaned)
  matched <- matched[!is.na(matched)]
  if (length(matched) == 0) {
    return(NA_integer_)
  }
  matched[1]
}

split_gene_values <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values)]
  pieces <- unlist(strsplit(values, "[,;|]", perl = TRUE), use.names = FALSE)
  trimws(pieces)
}

read_gmt_gene_sets <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("The GMT gene set file is empty.", call. = FALSE)
  }

  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) {
      return(NULL)
    }
    term_id <- trimws(parts[1])
    term_name <- trimws(parts[2])
    genes <- trimws(parts[-c(1, 2)])
    genes <- genes[nzchar(genes)]
    if (!nzchar(term_id) || length(genes) == 0) {
      return(NULL)
    }
    data.frame(
      term_id = term_id,
      term_name = if (nzchar(term_name)) term_name else term_id,
      gene_id = genes,
      stringsAsFactors = FALSE
    )
  })

  gene_sets <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    stop("No valid gene sets were found in the GMT file.", call. = FALSE)
  }
  unique(gene_sets)
}

read_tabular_gene_sets <- function(path, name = basename(path)) {
  data <- read_dge_table(path, name)
  term_col <- find_named_column(data, c("termid", "term", "geneset", "genesetid", "pathway", "pathwayid", "category", "name"))
  gene_col <- find_named_column(data, c("geneid", "gene", "genesymbol", "symbol", "featureid", "feature", "ensembl", "id"))
  term_name_col <- find_named_column(data, c("termname", "description", "genesetname", "pathwayname", "title"))

  if (!is.na(term_col) && !is.na(gene_col)) {
    rows <- lapply(seq_len(nrow(data)), function(row_index) {
      term_id <- trimws(as.character(data[[term_col]][row_index]))
      genes <- split_gene_values(data[[gene_col]][row_index])
      genes <- genes[nzchar(genes)]
      if (is.na(term_id) || !nzchar(term_id) || length(genes) == 0) {
        return(NULL)
      }
      term_name <- if (!is.na(term_name_col)) trimws(as.character(data[[term_name_col]][row_index])) else term_id
      data.frame(
        term_id = term_id,
        term_name = if (!is.na(term_name) && nzchar(term_name)) term_name else term_id,
        gene_id = genes,
        stringsAsFactors = FALSE
      )
    })
    gene_sets <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  } else {
    term_col <- if (is.na(term_col)) 1 else term_col
    term_name_col <- if (ncol(data) >= 3) 2 else NA_integer_
    gene_cols <- if (ncol(data) >= 3) setdiff(seq_len(ncol(data)), c(term_col, term_name_col)) else setdiff(seq_len(ncol(data)), term_col)

    rows <- lapply(seq_len(nrow(data)), function(row_index) {
      term_id <- trimws(as.character(data[[term_col]][row_index]))
      term_name <- if (!is.na(term_name_col)) trimws(as.character(data[[term_name_col]][row_index])) else term_id
      genes <- split_gene_values(unlist(data[row_index, gene_cols, drop = TRUE], use.names = FALSE))
      genes <- genes[nzchar(genes)]
      if (is.na(term_id) || !nzchar(term_id) || length(genes) == 0) {
        return(NULL)
      }
      data.frame(
        term_id = term_id,
        term_name = if (!is.na(term_name) && nzchar(term_name)) term_name else term_id,
        gene_id = genes,
        stringsAsFactors = FALSE
      )
    })
    gene_sets <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  }

  if (is.null(gene_sets) || nrow(gene_sets) == 0) {
    stop("No valid gene set mappings were found. Use GMT format or a table with term and gene columns.", call. = FALSE)
  }

  gene_sets <- unique(gene_sets)
  rownames(gene_sets) <- NULL
  gene_sets
}

read_gene_set_file <- function(path, name = basename(path)) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Gene set file could not be found.", call. = FALSE)
  }

  if (grepl("\\.gmt$", tolower(name))) {
    return(read_gmt_gene_sets(path))
  }

  read_tabular_gene_sets(path, name)
}

normalize_gene_ids <- function(gene_ids, case_sensitive = FALSE) {
  normalized <- trimws(as.character(gene_ids))
  normalized <- normalized[!is.na(normalized) & normalized != ""]
  if (!case_sensitive) {
    normalized <- toupper(normalized)
  }
  normalized
}

select_de_genes <- function(result_table, alpha = 0.05, lfc_cutoff = 1, gene_list = "both") {
  if (!all(c("gene_id", "log2FoldChange") %in% names(result_table))) {
    stop("The result table must contain gene_id and log2FoldChange columns.", call. = FALSE)
  }

  has_padj <- "padj" %in% names(result_table) && any(!is.na(result_table$padj))
  significant <- if (has_padj) {
    !is.na(result_table$padj) & result_table$padj < alpha
  } else {
    rep(TRUE, nrow(result_table))
  }

  lfc <- result_table$log2FoldChange
  up <- significant & !is.na(lfc) & lfc > lfc_cutoff
  down <- significant & !is.na(lfc) & lfc < -lfc_cutoff

  keep <- switch(
    gene_list,
    up = up,
    down = down,
    both = up | down,
    up | down
  )

  selected <- result_table[keep, , drop = FALSE]
  selected$regulation <- ifelse(selected$log2FoldChange > 0, "Up", "Down")
  selected <- selected[order(selected$regulation, selected$padj, -abs(selected$log2FoldChange), na.last = TRUE), , drop = FALSE]
  rownames(selected) <- NULL

  list(
    genes = unique(selected$gene_id),
    table = selected,
    has_padj = has_padj,
    selected_count = nrow(selected),
    up_count = sum(up),
    down_count = sum(down)
  )
}

run_ora_analysis <- function(
  universe_genes,
  selected_genes,
  gene_sets,
  min_set_size = 3,
  max_set_size = 2000,
  padj_cutoff = 0.1,
  case_sensitive = FALSE
) {
  if (min_set_size < 1) {
    stop("Minimum gene set size must be at least 1.", call. = FALSE)
  }
  if (max_set_size < min_set_size) {
    stop("Maximum gene set size must be greater than or equal to the minimum.", call. = FALSE)
  }
  if (!all(c("term_id", "term_name", "gene_id") %in% names(gene_sets))) {
    stop("Gene sets must contain term_id, term_name, and gene_id columns.", call. = FALSE)
  }

  universe_original <- trimws(as.character(universe_genes))
  universe_original <- universe_original[!is.na(universe_original) & universe_original != ""]
  universe_normalized <- normalize_gene_ids(universe_original, case_sensitive)
  keep_unique <- !duplicated(universe_normalized)
  universe_original <- universe_original[keep_unique]
  universe_normalized <- universe_normalized[keep_unique]

  selected_normalized <- normalize_gene_ids(selected_genes, case_sensitive)
  selected_normalized <- intersect(unique(selected_normalized), universe_normalized)

  if (length(universe_normalized) < 2) {
    stop("The ORA universe needs at least two genes.", call. = FALSE)
  }
  if (length(selected_normalized) == 0) {
    stop("No selected DEG genes are available for enrichment with the current thresholds.", call. = FALSE)
  }

  split_sets <- split(gene_sets, gene_sets$term_id)
  rows <- lapply(split_sets, function(term_data) {
    term_genes <- normalize_gene_ids(term_data$gene_id, case_sensitive)
    term_in_universe <- intersect(unique(term_genes), universe_normalized)
    set_size <- length(term_in_universe)
    if (set_size < min_set_size || set_size > max_set_size) {
      return(NULL)
    }

    overlap_normalized <- intersect(term_in_universe, selected_normalized)
    overlap_count <- length(overlap_normalized)
    universe_size <- length(universe_normalized)
    selected_count <- length(selected_normalized)
    expected <- selected_count * set_size / universe_size
    pvalue <- stats::phyper(
      q = overlap_count - 1,
      m = set_size,
      n = universe_size - set_size,
      k = selected_count,
      lower.tail = FALSE
    )
    overlap_genes <- universe_original[match(overlap_normalized, universe_normalized)]
    term_name <- term_data$term_name[which(nzchar(term_data$term_name))[1]]
    if (is.na(term_name) || !nzchar(term_name)) {
      term_name <- term_data$term_id[1]
    }

    data.frame(
      term_id = term_data$term_id[1],
      term_name = term_name,
      set_size = set_size,
      selected_genes = selected_count,
      overlap = overlap_count,
      expected = expected,
      fold_enrichment = if (expected > 0) overlap_count / expected else NA_real_,
      pvalue = pvalue,
      genes = paste(overlap_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  ora_table <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(ora_table) || nrow(ora_table) == 0) {
    stop("No gene sets overlapped the analysis universe after size filtering. Check gene IDs or change the size filters.", call. = FALSE)
  }

  ora_table$padj <- stats::p.adjust(ora_table$pvalue, method = "BH")
  ora_table <- ora_table[order(ora_table$padj, ora_table$pvalue, -ora_table$overlap), , drop = FALSE]
  ora_table <- ora_table[c(
    "term_id",
    "term_name",
    "set_size",
    "selected_genes",
    "overlap",
    "expected",
    "fold_enrichment",
    "pvalue",
    "padj",
    "genes"
  )]
  rownames(ora_table) <- NULL

  summary_table <- data.frame(
    metric = c("Universe genes", "Selected genes", "Gene sets loaded", "Gene sets tested", "Significant terms"),
    value = c(
      length(universe_normalized),
      length(selected_normalized),
      length(split_sets),
      nrow(ora_table),
      sum(!is.na(ora_table$padj) & ora_table$padj < padj_cutoff)
    ),
    stringsAsFactors = FALSE
  )

  list(
    results = ora_table,
    summary = summary_table,
    selected_genes = selected_normalized,
    universe_genes = universe_normalized,
    min_set_size = min_set_size,
    max_set_size = max_set_size,
    padj_cutoff = padj_cutoff,
    case_sensitive = case_sensitive
  )
}

strip_kegg_species_suffix <- function(term_names) {
  sub("\\s+-\\s+[^-]+\\s*\\([^)]*\\)\\s*$", "", term_names)
}

annotation_gene_ids_by_biotype <- function(result_table, annotation_info = NULL, match_mode = "auto", biotype = "protein_coding") {
  if (is.null(annotation_info)) {
    result_ids <- unique(trimws(as.character(result_table$gene_id)))
    result_ids <- result_ids[!is.na(result_ids) & nzchar(result_ids)]
    return(list(
      ids = result_ids,
      universe_ids = result_ids,
      matched = length(result_ids),
      total = length(result_ids),
      match_mode = "none"
    ))
  }

  annotation <- annotation_info$table
  resolved_mode <- resolve_annotation_match_mode(result_table, annotation, match_mode)
  annotation_key_col <- if (identical(resolved_mode, "symbol")) "gene_symbol" else "ensembl_gene_id"
  keep <- !is.na(annotation$gene_biotype) &
    annotation$gene_biotype == biotype &
    !is.na(annotation[[annotation_key_col]]) &
    annotation[[annotation_key_col]] != ""

  allowed <- normalize_gene_ids(annotation[[annotation_key_col]][keep], case_sensitive = FALSE)
  allowed <- unique(allowed[nzchar(allowed)])
  result_ids <- trimws(as.character(result_table$gene_id))
  result_normalized <- normalize_gene_ids(result_ids, case_sensitive = FALSE)
  matched <- result_ids[result_normalized %in% allowed]

  list(
    ids = unique(matched[!is.na(matched) & matched != ""]),
    universe_ids = allowed,
    matched = length(unique(matched[!is.na(matched) & matched != ""])),
    total = length(allowed),
    match_mode = resolved_mode
  )
}

run_rosbys_lab_enrichment <- function(
  result_table,
  gene_sets,
  annotation_info = NULL,
  alpha = 0.1,
  lfc_cutoff = log2(2),
  padj_cutoff = 0.01,
  max_terms_per_direction = 5,
  min_overlap = 2,
  case_sensitive = FALSE
) {
  if (!all(c("gene_id", "log2FoldChange") %in% names(result_table))) {
    stop("The result table must contain gene_id and log2FoldChange columns.", call. = FALSE)
  }
  if (!all(c("term_id", "term_name", "gene_id") %in% names(gene_sets))) {
    stop("Gene sets must contain term_id, term_name, and gene_id columns.", call. = FALSE)
  }

  has_padj <- "padj" %in% names(result_table) && any(!is.na(result_table$padj))
  significant <- if (has_padj) {
    !is.na(result_table$padj) & result_table$padj <= alpha
  } else {
    rep(TRUE, nrow(result_table))
  }

  lfc <- result_table$log2FoldChange
  protein_coding <- annotation_gene_ids_by_biotype(result_table, annotation_info, "auto", "protein_coding")
  protein_coding_norm <- normalize_gene_ids(protein_coding$ids, case_sensitive)
  universe_norm <- unique(normalize_gene_ids(protein_coding$universe_ids, case_sensitive))
  universe_norm <- universe_norm[!is.na(universe_norm) & nzchar(universe_norm)]
  result_norm <- normalize_gene_ids(result_table$gene_id, case_sensitive)
  protein_coding_result <- result_norm %in% protein_coding_norm

  direction_specs <- list(
    list(group = "Upregulated", keep = significant & protein_coding_result & !is.na(lfc) & lfc >= lfc_cutoff),
    list(group = "Downregulated", keep = significant & protein_coding_result & !is.na(lfc) & lfc <= -lfc_cutoff)
  )

  split_sets <- split(gene_sets, gene_sets$term_id)
  universe_size <- length(universe_norm)
  if (universe_size < 2) {
    stop("The Rosby's Lab-style enrichment universe needs at least two genes.", call. = FALSE)
  }

  all_rows <- list()
  selected_tables <- list()
  tested_counts <- integer()

  for (direction in direction_specs) {
    selected <- result_table[direction$keep, , drop = FALSE]
    selected$regulation <- direction$group
    selected_tables[[direction$group]] <- selected
    selected_norm <- unique(normalize_gene_ids(selected$gene_id, case_sensitive))
    selected_norm <- selected_norm[nzchar(selected_norm)]
    selected_count <- length(selected_norm)
    selected_lookup <- data.frame(
      normalized = normalize_gene_ids(selected$gene_id, case_sensitive),
      gene_id = trimws(as.character(selected$gene_id)),
      stringsAsFactors = FALSE
    )
    selected_lookup <- selected_lookup[
      !is.na(selected_lookup$normalized) &
        nzchar(selected_lookup$normalized) &
        !is.na(selected_lookup$gene_id) &
        nzchar(selected_lookup$gene_id),
      ,
      drop = FALSE
    ]
    selected_lookup <- selected_lookup[!duplicated(selected_lookup$normalized), , drop = FALSE]

    if (selected_count <= 2) {
      tested_counts[[direction$group]] <- 0L
      next
    }

    rows <- lapply(split_sets, function(term_data) {
      term_genes <- unique(normalize_gene_ids(term_data$gene_id, case_sensitive))
      term_genes <- term_genes[!is.na(term_genes) & nzchar(term_genes)]
      term_genes <- intersect(term_genes, universe_norm)
      set_size <- length(term_genes)
      overlap_normalized <- intersect(term_genes, selected_norm)
      overlap_count <- length(overlap_normalized)
      if (overlap_count < min_overlap || set_size < min_overlap || set_size > universe_size) {
        return(NULL)
      }

      pvalue <- stats::phyper(
        q = overlap_count - 1,
        m = selected_count,
        n = universe_size - selected_count,
        k = set_size,
        lower.tail = FALSE
      )
      term_name <- term_data$term_name[which(nzchar(term_data$term_name))[1]]
      if (is.na(term_name) || !nzchar(term_name)) {
        term_name <- term_data$term_id[1]
      }

      overlap_genes <- selected_lookup$gene_id[match(overlap_normalized, selected_lookup$normalized)]
      overlap_genes <- overlap_genes[!is.na(overlap_genes) & overlap_genes != ""]
      data.frame(
        group = direction$group,
        term_id = term_data$term_id[1],
        term_name = term_name,
        set_size = set_size,
        selected_genes = selected_count,
        overlap = overlap_count,
        expected = selected_count * set_size / universe_size,
        fold_enrichment = overlap_count / (selected_count * set_size / universe_size),
        pvalue = pvalue,
        genes = paste(unique(overlap_genes), collapse = ";"),
        stringsAsFactors = FALSE
      )
    })

    rows <- rows[!vapply(rows, is.null, logical(1))]
    direction_table <- if (length(rows) > 0) do.call(rbind, rows) else NULL
    if (is.null(direction_table) || nrow(direction_table) == 0) {
      tested_counts[[direction$group]] <- 0L
      next
    }

    direction_table$padj <- stats::p.adjust(direction_table$pvalue, method = "fdr")
    direction_table <- direction_table[order(direction_table$padj, direction_table$pvalue), , drop = FALSE]
    tested_counts[[direction$group]] <- nrow(direction_table)
    direction_table <- direction_table[!is.na(direction_table$padj) & direction_table$padj < padj_cutoff, , drop = FALSE]
    if (nrow(direction_table) > max_terms_per_direction) {
      direction_table <- direction_table[seq_len(max_terms_per_direction), , drop = FALSE]
    }
    all_rows[[direction$group]] <- direction_table
  }

  all_rows <- all_rows[!vapply(all_rows, is.null, logical(1))]
  enrichment_table <- if (length(all_rows) > 0) do.call(rbind, all_rows) else NULL
  if (is.null(enrichment_table) || nrow(enrichment_table) == 0) {
    stop("No significant Rosby's Lab-style enrichment found with the current DEG and FDR thresholds.", call. = FALSE)
  }
  rownames(enrichment_table) <- NULL
  enrichment_table <- enrichment_table[c(
    "group",
    "term_id",
    "term_name",
    "set_size",
    "selected_genes",
    "overlap",
    "expected",
    "fold_enrichment",
    "pvalue",
    "padj",
    "genes"
  )]

  selected_table <- do.call(rbind, selected_tables[!vapply(selected_tables, is.null, logical(1))])
  if (is.null(selected_table)) {
    selected_table <- result_table[FALSE, , drop = FALSE]
  }
  rownames(selected_table) <- NULL

  summary_table <- data.frame(
    metric = c(
      "Mode",
      "Universe genes",
      "Upregulated DEG genes",
      "Downregulated DEG genes",
      "Gene sets loaded",
      "Gene sets tested",
      "Displayed terms",
      "Gene-set source"
    ),
    value = c(
      "Rosby's Lab-style ORA",
      universe_size,
      nrow(selected_tables[["Upregulated"]]),
      nrow(selected_tables[["Downregulated"]]),
      length(split_sets),
      sum(tested_counts),
      nrow(enrichment_table),
      "Selected GO, KEGG, or custom source"
    ),
    stringsAsFactors = FALSE
  )

  list(
    mode = "rosbys_lab",
    results = enrichment_table,
    summary = summary_table,
    selected_table = selected_table,
    selected_genes = normalize_gene_ids(selected_table$gene_id, case_sensitive),
    universe_genes = universe_norm,
    padj_cutoff = padj_cutoff,
    max_terms_per_direction = max_terms_per_direction,
    min_overlap = min_overlap,
    has_padj = has_padj,
    protein_coding_summary = protein_coding
  )
}

format_enrichment_preview <- function(ora_table, n = 100) {
  preview <- utils::head(ora_table, n)
  pvalue_cols <- intersect(c("pvalue", "padj"), names(preview))
  for (column_name in pvalue_cols) {
    preview[[column_name]] <- ifelse(
      is.na(preview[[column_name]]),
      "NA",
      formatC(preview[[column_name]], format = "e", digits = 3)
    )
  }
  numeric_cols <- intersect(c("expected", "fold_enrichment"), names(preview))
  preview[numeric_cols] <- lapply(preview[numeric_cols], function(column) signif(column, 4))
  preview
}

format_rosbys_lab_enrichment_export <- function(ora_table) {
  term_url <- ifelse(
    grepl("^[a-z][a-z][a-z][0-9]{5}$", ora_table$term_id),
    paste0("http://www.genome.jp/kegg-bin/show_pathway?", ora_table$term_id, " "),
    ""
  )
  data.frame(
    group = ora_table$group,
    FDR = ora_table$padj,
    nGenes = ora_table$overlap,
    "Pathway size" = ora_table$set_size,
    "Fold enriched" = ora_table$fold_enrichment,
    Pathway = strip_kegg_species_suffix(ora_table$term_name),
    URL = term_url,
    Genes = gsub(";", ", ", ora_table$genes, fixed = TRUE),
    check.names = FALSE
  )
}

format_rosbys_lab_enrichment_preview <- function(ora_table, n = 100) {
  preview <- utils::head(ora_table, n)
  data.frame(
    "Grp." = ifelse(duplicated(preview$group), "", preview$group),
    "Adj.Pval" = formatC(preview$padj, format = "e", digits = 2),
    "nGenes" = preview$overlap,
    "Pathway size" = preview$set_size,
    "Fold" = round(preview$fold_enrichment, 1),
    "Pathway" = strip_kegg_species_suffix(preview$term_name),
    "Genes" = preview$genes,
    check.names = FALSE
  )
}

make_ora_plot <- function(ora_table, top_n = 20, padj_cutoff = 0.1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }
  validate_rows <- ora_table[ora_table$overlap > 0, , drop = FALSE]
  if (nrow(validate_rows) == 0) {
    stop("No gene sets have overlapping selected genes.", call. = FALSE)
  }

  plot_data <- utils::head(validate_rows[order(validate_rows$padj, validate_rows$pvalue), , drop = FALSE], top_n)
  plot_data$score <- -log10(pmax(plot_data$padj, .Machine$double.xmin))
  plot_data$term_label <- ifelse(nzchar(plot_data$term_name), plot_data$term_name, plot_data$term_id)
  plot_data$term_label <- factor(plot_data$term_label, levels = rev(plot_data$term_label))
  plot_data$significant <- ifelse(plot_data$padj < padj_cutoff, "FDR pass", "FDR above cutoff")

  ggplot2::ggplot(plot_data, ggplot2::aes(x = score, y = term_label, fill = significant)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(
      ggplot2::aes(label = overlap),
      hjust = -0.2,
      size = 3.3
    ) +
    ggplot2::scale_fill_manual(values = c("FDR pass" = "#0F766E", "FDR above cutoff" = "#6B7280")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      x = "-log10 adjusted p-value",
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

prepare_pathway_ranks <- function(
  result_table,
  rank_metric = "log2fc",
  gene_padj_cutoff = 1,
  absolute_ranking = FALSE,
  case_sensitive = FALSE
) {
  if (!all(c("gene_id", "log2FoldChange") %in% names(result_table))) {
    stop("Pathway analysis needs gene_id and log2FoldChange columns.", call. = FALSE)
  }
  if (!rank_metric %in% c("log2fc", "signed_p")) {
    stop("Choose a valid pathway ranking metric.", call. = FALSE)
  }
  if (!is.numeric(gene_padj_cutoff) || length(gene_padj_cutoff) != 1 ||
      is.na(gene_padj_cutoff) || gene_padj_cutoff <= 0 || gene_padj_cutoff > 1) {
    stop("The gene adjusted p-value filter must be greater than 0 and at most 1.", call. = FALSE)
  }

  gene_ids <- trimws(as.character(result_table$gene_id))
  normalized_ids <- gene_ids
  if (!case_sensitive) {
    normalized_ids <- toupper(normalized_ids)
  }
  log2fc <- suppressWarnings(as.numeric(result_table$log2FoldChange))
  pvalue <- if ("pvalue" %in% names(result_table)) {
    suppressWarnings(as.numeric(result_table$pvalue))
  } else {
    rep(NA_real_, nrow(result_table))
  }
  padj <- if ("padj" %in% names(result_table)) {
    suppressWarnings(as.numeric(result_table$padj))
  } else {
    rep(NA_real_, nrow(result_table))
  }

  if (identical(rank_metric, "signed_p")) {
    if (!any(is.finite(pvalue))) {
      stop("Signed p-value ranking needs a result table with p-values.", call. = FALSE)
    }
    rank_score <- sign(log2fc) * -log10(pmax(pvalue, .Machine$double.xmin))
  } else {
    rank_score <- log2fc
  }

  keep <- !is.na(gene_ids) & gene_ids != "" &
    !is.na(normalized_ids) & normalized_ids != "" &
    is.finite(rank_score)
  if (gene_padj_cutoff < 1) {
    if (!any(is.finite(padj))) {
      stop("The gene FDR filter cannot be used because adjusted p-values are unavailable.", call. = FALSE)
    }
    keep <- keep & is.finite(padj) & padj <= gene_padj_cutoff
  }

  ranked_table <- data.frame(
    gene_id = gene_ids[keep],
    normalized_gene_id = normalized_ids[keep],
    rank_score = rank_score[keep],
    log2FoldChange = log2fc[keep],
    pvalue = pvalue[keep],
    padj = padj[keep],
    stringsAsFactors = FALSE
  )
  if (nrow(ranked_table) < 10) {
    stop("Pathway analysis needs at least 10 ranked genes after filtering.", call. = FALSE)
  }

  ranked_table <- ranked_table[
    order(abs(ranked_table$rank_score), decreasing = TRUE, na.last = NA),
    ,
    drop = FALSE
  ]
  ranked_table <- ranked_table[!duplicated(ranked_table$normalized_gene_id), , drop = FALSE]
  if (isTRUE(absolute_ranking)) {
    ranked_table$rank_score <- abs(ranked_table$rank_score)
  }
  ranked_table <- ranked_table[
    order(ranked_table$rank_score, decreasing = TRUE, na.last = NA),
    ,
    drop = FALSE
  ]
  rownames(ranked_table) <- NULL

  ranked_stats <- ranked_table$rank_score
  names(ranked_stats) <- ranked_table$normalized_gene_id
  list(table = ranked_table, stats = ranked_stats)
}

prepare_pathway_collections <- function(gene_sets, ranked_gene_ids, min_set_size = 10, max_set_size = 500, case_sensitive = FALSE) {
  if (!all(c("term_id", "term_name", "gene_id") %in% names(gene_sets))) {
    stop("Gene sets must contain term_id, term_name, and gene_id columns.", call. = FALSE)
  }
  if (min_set_size < 1) {
    stop("Minimum pathway size must be at least 1.", call. = FALSE)
  }
  if (max_set_size < min_set_size) {
    stop("Maximum pathway size must be greater than or equal to the minimum.", call. = FALSE)
  }

  mapped_sets <- gene_sets
  mapped_sets$gene_id <- trimws(as.character(mapped_sets$gene_id))
  if (!case_sensitive) {
    mapped_sets$gene_id <- toupper(mapped_sets$gene_id)
  }
  mapped_sets <- mapped_sets[
    !is.na(mapped_sets$term_id) & mapped_sets$term_id != "" &
      !is.na(mapped_sets$gene_id) & mapped_sets$gene_id != "",
    ,
    drop = FALSE
  ]

  split_sets <- split(mapped_sets, mapped_sets$term_id)
  pathway_list <- lapply(split_sets, function(term_data) {
    intersect(unique(term_data$gene_id), ranked_gene_ids)
  })
  sizes <- vapply(pathway_list, length, integer(1))
  pathway_list <- pathway_list[sizes >= min_set_size & sizes <= max_set_size]
  if (length(pathway_list) == 0) {
    stop("No pathways overlap enough ranked genes after size filtering. Check the organism, gene IDs, or pathway-size limits.", call. = FALSE)
  }

  term_rows <- lapply(names(pathway_list), function(term_id) {
    term_data <- split_sets[[term_id]]
    term_name <- term_data$term_name[which(!is.na(term_data$term_name) & nzchar(term_data$term_name))[1]]
    if (is.na(term_name) || !nzchar(term_name)) {
      term_name <- term_id
    }
    data.frame(
      term_id = term_id,
      term_name = term_name,
      stringsAsFactors = FALSE
    )
  })

  list(
    pathways = pathway_list,
    terms = do.call(rbind, term_rows),
    loaded_count = length(split_sets)
  )
}

run_preranked_pathway_analysis <- function(
  result_table,
  gene_sets,
  rank_metric = "log2fc",
  min_set_size = 10,
  max_set_size = 500,
  padj_cutoff = 0.05,
  gene_padj_cutoff = 1,
  absolute_ranking = FALSE,
  case_sensitive = FALSE
) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("The fgsea package is required for Pathway Analysis. Run scripts/install_packages.R, then restart the app.", call. = FALSE)
  }
  if (!is.numeric(padj_cutoff) || length(padj_cutoff) != 1 ||
      is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
    stop("The pathway FDR cutoff must be greater than 0 and at most 1.", call. = FALSE)
  }

  ranks <- prepare_pathway_ranks(
    result_table = result_table,
    rank_metric = rank_metric,
    gene_padj_cutoff = gene_padj_cutoff,
    absolute_ranking = absolute_ranking,
    case_sensitive = case_sensitive
  )
  collections <- prepare_pathway_collections(
    gene_sets = gene_sets,
    ranked_gene_ids = names(ranks$stats),
    min_set_size = min_set_size,
    max_set_size = max_set_size,
    case_sensitive = case_sensitive
  )
  score_type <- if (isTRUE(absolute_ranking) || all(ranks$stats >= 0)) {
    "pos"
  } else if (all(ranks$stats <= 0)) {
    "neg"
  } else {
    "std"
  }

  set.seed(42)
  fgsea_result <- suppressWarnings(
    fgsea::fgseaMultilevel(
      pathways = collections$pathways,
      stats = ranks$stats,
      minSize = min_set_size,
      maxSize = max_set_size,
      eps = 0,
      scoreType = score_type,
      nproc = 1
    )
  )
  fgsea_result <- as.data.frame(fgsea_result, stringsAsFactors = FALSE)
  if (nrow(fgsea_result) == 0) {
    stop("No pathways could be tested with the current ranked genes and size filters.", call. = FALSE)
  }

  term_match <- match(fgsea_result$pathway, collections$terms$term_id)
  normalized_to_original <- stats::setNames(ranks$table$gene_id, ranks$table$normalized_gene_id)
  leading_edges <- lapply(fgsea_result$leadingEdge, function(ids) {
    original <- unname(normalized_to_original[as.character(ids)])
    original <- original[!is.na(original) & original != ""]
    unique(original)
  })
  pathway_table <- data.frame(
    term_id = fgsea_result$pathway,
    term_name = collections$terms$term_name[term_match],
    direction = ifelse(fgsea_result$NES >= 0, "Positive", "Negative"),
    set_size = fgsea_result$size,
    ES = fgsea_result$ES,
    NES = fgsea_result$NES,
    pvalue = fgsea_result$pval,
    padj = fgsea_result$padj,
    leading_edge_count = vapply(leading_edges, length, integer(1)),
    leading_edge_genes = vapply(leading_edges, paste, collapse = ";", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  )
  pathway_table$term_name[
    is.na(pathway_table$term_name) | pathway_table$term_name == ""
  ] <- pathway_table$term_id[
    is.na(pathway_table$term_name) | pathway_table$term_name == ""
  ]
  pathway_table <- pathway_table[
    order(pathway_table$padj, pathway_table$pvalue, -abs(pathway_table$NES), na.last = TRUE),
    ,
    drop = FALSE
  ]
  rownames(pathway_table) <- NULL

  summary_table <- data.frame(
    metric = c("Ranked genes", "Gene sets loaded", "Pathways tested", "Significant pathways", "Positive NES", "Negative NES"),
    value = c(
      length(ranks$stats),
      collections$loaded_count,
      nrow(pathway_table),
      sum(!is.na(pathway_table$padj) & pathway_table$padj < padj_cutoff),
      sum(pathway_table$NES > 0, na.rm = TRUE),
      sum(pathway_table$NES < 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  list(
    mode = "preranked_gsea",
    results = pathway_table,
    summary = summary_table,
    ranked_table = ranks$table,
    ranked_stats = ranks$stats,
    pathways = collections$pathways,
    rank_metric = rank_metric,
    min_set_size = min_set_size,
    max_set_size = max_set_size,
    padj_cutoff = padj_cutoff,
    gene_padj_cutoff = gene_padj_cutoff,
    absolute_ranking = absolute_ranking,
    score_type = score_type
  )
}

format_pathway_preview <- function(pathway_table, n = 100, show_ids = FALSE) {
  preview <- utils::head(pathway_table, n)
  if (!isTRUE(show_ids)) {
    preview$term_id <- NULL
  }
  for (column_name in intersect(c("pvalue", "padj"), names(preview))) {
    preview[[column_name]] <- ifelse(
      is.na(preview[[column_name]]),
      "NA",
      formatC(preview[[column_name]], format = "e", digits = 3)
    )
  }
  for (column_name in intersect(c("ES", "NES"), names(preview))) {
    preview[[column_name]] <- round(preview[[column_name]], 3)
  }
  display_names <- c(
    term_id = "Pathway ID",
    term_name = "Pathway",
    direction = "Direction",
    set_size = "Set size",
    ES = "ES",
    NES = "NES",
    pvalue = "P-value",
    padj = "Adjusted p-value",
    leading_edge_count = "Leading-edge count",
    leading_edge_genes = "Leading-edge genes"
  )
  names(preview) <- ifelse(
    names(preview) %in% names(display_names),
    unname(display_names[names(preview)]),
    names(preview)
  )
  preview
}

make_pathway_summary_plot <- function(pathway_table, top_n = 20, padj_cutoff = 0.05, show_ids = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }
  plot_data <- pathway_table[is.finite(pathway_table$NES), , drop = FALSE]
  if (nrow(plot_data) == 0) {
    stop("No pathway enrichment scores are available to plot.", call. = FALSE)
  }
  significant <- plot_data[!is.na(plot_data$padj) & plot_data$padj < padj_cutoff, , drop = FALSE]
  if (nrow(significant) > 0) {
    plot_data <- significant
  }
  plot_data <- utils::head(
    plot_data[order(plot_data$padj, plot_data$pvalue, -abs(plot_data$NES), na.last = TRUE), , drop = FALSE],
    top_n
  )
  labels <- if (isTRUE(show_ids)) {
    paste0(plot_data$term_name, " [", plot_data$term_id, "]")
  } else {
    plot_data$term_name
  }
  labels <- ifelse(nchar(labels) > 70, paste0(substr(labels, 1, 67), "..."), labels)
  plot_data$term_label <- factor(make.unique(labels), levels = rev(make.unique(labels)))
  plot_data$fdr_score <- -log10(pmax(plot_data$padj, .Machine$double.xmin))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = NES, y = term_label, fill = direction)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, color = "#64748B", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("Positive" = "#B42318", "Negative" = "#2563EB")) +
    ggplot2::labs(
      x = "Normalized enrichment score (NES)",
      y = NULL,
      fill = "Enrichment"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

make_pathway_enrichment_plot <- function(pathway_result, term_id) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("The fgsea package is required for pathway enrichment plots.", call. = FALSE)
  }
  if (is.null(term_id) || !term_id %in% names(pathway_result$pathways)) {
    stop("Choose a pathway to display its running enrichment score.", call. = FALSE)
  }
  result_row <- pathway_result$results[pathway_result$results$term_id == term_id, , drop = FALSE]
  title <- if (nrow(result_row) > 0) result_row$term_name[1] else term_id
  fgsea::plotEnrichment(pathway_result$pathways[[term_id]], pathway_result$ranked_stats) +
    ggplot2::labs(
      title = title,
      subtitle = if (nrow(result_row) > 0) {
        sprintf("NES %.2f | adjusted p %s", result_row$NES[1], formatC(result_row$padj[1], format = "e", digits = 2))
      } else {
        NULL
      },
      x = "Ranked gene position",
      y = "Running enrichment score"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

pathway_leading_edge_table <- function(pathway_result, term_id) {
  row <- pathway_result$results[pathway_result$results$term_id == term_id, , drop = FALSE]
  if (nrow(row) == 0 || !nzchar(row$leading_edge_genes[1])) {
    return(pathway_result$ranked_table[FALSE, , drop = FALSE])
  }
  genes <- strsplit(row$leading_edge_genes[1], ";", fixed = TRUE)[[1]]
  match_index <- match(genes, pathway_result$ranked_table$gene_id)
  leading <- pathway_result$ranked_table[match_index[!is.na(match_index)], , drop = FALSE]
  leading <- leading[c("gene_id", "rank_score", "log2FoldChange", "pvalue", "padj")]
  rownames(leading) <- NULL
  leading
}

make_pathway_heatmap <- function(workflow_result, gene_ids, max_genes = 40) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }
  expression_table <- workflow_result$normalized_counts
  if (is.null(expression_table) || !"gene_id" %in% names(expression_table)) {
    stop("A pathway heatmap is available only when sample-level expression values were analyzed.", call. = FALSE)
  }

  expression_ids <- toupper(trimws(as.character(expression_table$gene_id)))
  requested_ids <- toupper(trimws(as.character(gene_ids)))
  matched <- match(requested_ids, expression_ids)
  matched <- matched[!is.na(matched)]
  matched <- unique(matched)
  if (length(matched) == 0) {
    stop("No leading-edge genes matched the sample-level expression matrix.", call. = FALSE)
  }
  matched <- utils::head(matched, max_genes)
  values <- as.matrix(expression_table[matched, setdiff(names(expression_table), "gene_id"), drop = FALSE])
  storage.mode(values) <- "double"
  if (identical(workflow_result$input_type, "counts")) {
    values <- log2(values + 1)
  }
  rownames(values) <- expression_table$gene_id[matched]
  row_sd <- apply(values, 1, stats::sd, na.rm = TRUE)
  values <- values[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
  if (nrow(values) == 0) {
    stop("Leading-edge genes have no variable sample-level expression to display.", call. = FALSE)
  }
  scaled <- t(scale(t(values)))
  scaled[!is.finite(scaled)] <- 0
  heatmap_data <- as.data.frame(as.table(scaled), stringsAsFactors = FALSE)
  names(heatmap_data) <- c("gene_id", "sample_id", "z_score")
  heatmap_data$gene_id <- factor(heatmap_data$gene_id, levels = rev(rownames(scaled)))
  heatmap_data$sample_id <- factor(heatmap_data$sample_id, levels = colnames(scaled))

  ggplot2::ggplot(heatmap_data, ggplot2::aes(x = sample_id, y = gene_id, fill = z_score)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#2563EB",
      mid = "#F8FAFC",
      high = "#B42318",
      midpoint = 0,
      limits = c(-max(abs(heatmap_data$z_score)), max(abs(heatmap_data$z_score)))
    ) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Row z-score") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

standardize_fold_change_table <- function(data, require_adjusted_p = FALSE, fold_change_scale = "auto") {
  if (ncol(data) < 2) {
    stop("Fold-change input needs at least a gene ID column and a fold-change column.", call. = FALSE)
  }

  cleaned <- clean_column_name(names(data))
  gene_col <- find_named_column(data, c("geneid", "gene", "id", "genesymbol", "symbol", "featureid", "feature"))
  if (is.na(gene_col)) {
    gene_col <- 1
  }

  log_fc_col <- find_named_column(data, c("log2foldchange", "log2fc", "logfc", "lfc", "avglog2fc", "avglogfc"))
  linear_fc_col <- find_named_column(data, c("foldchange", "foldchangevalue", "fc"))

  fc_col <- if (!is.na(log_fc_col)) log_fc_col else linear_fc_col
  if (is.na(fc_col)) {
    numeric_cols <- which(vapply(data, function(column) {
      suppressWarnings(!all(is.na(as.numeric(as.character(column)))))
    }, logical(1)))
    numeric_cols <- setdiff(numeric_cols, gene_col)
    if (length(numeric_cols) == 0) {
      stop("Could not find a numeric fold-change column.", call. = FALSE)
    }
    fc_col <- numeric_cols[1]
  }

  padj_col <- find_named_column(data, c("padj", "adjustedpvalue", "adjpvalue", "pvalueadjusted", "qvalue", "fdr"))
  pvalue_col <- find_named_column(data, c("pvalue", "pval", "rawpvalue", "p"))
  base_mean_col <- find_named_column(data, c("basemean", "mean", "average", "averageexpression", "meanexpression"))

  if (require_adjusted_p && is.na(padj_col)) {
    stop("Could not find an adjusted p-value column. Accepted names include padj, FDR, qvalue, and adjusted_p_value.", call. = FALSE)
  }

  gene_ids <- trimws(as.character(data[[gene_col]]))
  if (any(is.na(gene_ids) | gene_ids == "")) {
    stop("The gene ID column contains empty values.", call. = FALSE)
  }

  raw_fc <- suppressWarnings(as.numeric(gsub(",", "", as.character(data[[fc_col]]), fixed = TRUE)))
  if (anyNA(raw_fc)) {
    stop("The fold-change column contains missing or non-numeric values.", call. = FALSE)
  }

  resolved_scale <- fold_change_scale
  if (identical(resolved_scale, "auto")) {
    resolved_scale <- if (!is.na(log_fc_col) || any(raw_fc <= 0)) "log2" else "linear"
  }

  if (identical(resolved_scale, "linear")) {
    if (any(raw_fc <= 0)) {
      stop("Linear fold changes must be greater than zero. Choose log2 fold-change scale.", call. = FALSE)
    }
    log2_fold_change <- log2(raw_fc)
  } else {
    log2_fold_change <- raw_fc
  }

  pvalue <- if (is.na(pvalue_col)) rep(NA_real_, length(gene_ids)) else suppressWarnings(as.numeric(data[[pvalue_col]]))
  padj <- if (is.na(padj_col)) rep(NA_real_, length(gene_ids)) else suppressWarnings(as.numeric(data[[padj_col]]))
  base_mean <- if (is.na(base_mean_col)) rep(NA_real_, length(gene_ids)) else suppressWarnings(as.numeric(data[[base_mean_col]]))

  result_table <- data.frame(
    gene_id = gene_ids,
    baseMean = base_mean,
    log2FoldChange = log2_fold_change,
    lfcSE = NA_real_,
    stat = NA_real_,
    pvalue = pvalue,
    padj = padj,
    stringsAsFactors = FALSE
  )

  list(
    result_table = sort_result_table(result_table),
    fold_change_column = names(data)[fc_col],
    adjusted_p_column = if (is.na(padj_col)) NULL else names(data)[padj_col],
    fold_change_scale = resolved_scale
  )
}

run_fold_change_workflow <- function(data, require_adjusted_p = FALSE, fold_change_scale = "auto", alpha = 0.05) {
  standardized <- standardize_fold_change_table(
    data = data,
    require_adjusted_p = require_adjusted_p,
    fold_change_scale = fold_change_scale
  )

  warnings <- c(
    "Fold-change mode does not run DESeq2. It visualizes and exports the fold-change table you supplied.",
    sprintf("Fold-change column used: %s.", standardized$fold_change_column),
    sprintf("Fold-change scale used: %s.", standardized$fold_change_scale)
  )
  if (is.null(standardized$adjusted_p_column)) {
    warnings <- c(warnings, "No adjusted p-value column was supplied, so significance calls are not available.")
  } else {
    warnings <- c(warnings, sprintf("Adjusted p-value column used: %s.", standardized$adjusted_p_column))
  }

  list(
    dds = NULL,
    deseq_result = NULL,
    results = standardized$result_table,
    normalized_counts = NULL,
    pca_data = NULL,
    percent_var = NULL,
    warnings = unique(warnings),
    genes_tested = nrow(standardized$result_table),
    samples_tested = NA_integer_,
    condition_col = NULL,
    treatment_level = NULL,
    reference_level = NULL,
    batch_col = NULL,
    alpha = alpha,
    input_type = if (require_adjusted_p) "fc_padj" else "fc_only",
    workflow_label = if (require_adjusted_p) "Fold changes and adjusted p-values" else "Fold changes only"
  )
}

make_pca_plot <- function(workflow_result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }

  ggplot2::ggplot(
    workflow_result$pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = condition, label = name)
  ) +
    ggplot2::geom_point(size = 3.2, alpha = 0.9) +
    ggplot2::geom_text(vjust = -0.8, size = 3.2, show.legend = FALSE) +
    ggplot2::labs(
      x = sprintf("PC1: %s%% variance", workflow_result$percent_var[1]),
      y = sprintf("PC2: %s%% variance", workflow_result$percent_var[2]),
      color = "Condition"
    ) +
    ggplot2::scale_color_manual(values = c("#0F766E", "#B42318", "#5B5BD6", "#8A5A00", "#406A38")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

make_volcano_plot <- function(result_table, alpha = 0.05, lfc_cutoff = 1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }

  plot_data <- result_table
  has_padj <- "padj" %in% names(plot_data) && any(!is.na(plot_data$padj))

  if (has_padj) {
    plot_data$neg_log10_padj <- -log10(pmax(plot_data$padj, .Machine$double.xmin, na.rm = FALSE))
    plot_data$status <- "Not significant"
    plot_data$status[!is.na(plot_data$padj) &
      plot_data$padj < alpha &
      plot_data$log2FoldChange > lfc_cutoff] <- "Up"
    plot_data$status[!is.na(plot_data$padj) &
      plot_data$padj < alpha &
      plot_data$log2FoldChange < -lfc_cutoff] <- "Down"

    return(
      ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = status)
      ) +
        ggplot2::geom_point(alpha = 0.75, size = 1.9, na.rm = TRUE) +
        ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "#6B7280") +
        ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "#6B7280") +
        ggplot2::labs(
          x = "log2 fold change",
          y = "-log10 adjusted p-value",
          color = NULL
        ) +
        ggplot2::scale_color_manual(
          values = c("Down" = "#2563EB", "Not significant" = "#6B7280", "Up" = "#B42318")
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "bottom"
        )
    )
  }

  plot_data$status <- "Below cutoff"
  plot_data$status[plot_data$log2FoldChange > lfc_cutoff] <- "Up"
  plot_data$status[plot_data$log2FoldChange < -lfc_cutoff] <- "Down"
  plot_data$absolute_lfc <- abs(plot_data$log2FoldChange)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = log2FoldChange, y = absolute_lfc, color = status)
  ) +
    ggplot2::geom_point(alpha = 0.75, size = 1.9, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "#6B7280") +
    ggplot2::labs(
      x = "log2 fold change",
      y = "|log2 fold change|",
      color = NULL
    ) +
    ggplot2::scale_color_manual(
      values = c("Down" = "#2563EB", "Below cutoff" = "#6B7280", "Up" = "#B42318")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

regulation_summary_counts <- function(result_table, alpha = 0.05, lfc_cutoff = 1) {
  has_padj <- "padj" %in% names(result_table) && any(!is.na(result_table$padj))
  significant <- if (has_padj) {
    !is.na(result_table$padj) & result_table$padj < alpha
  } else {
    rep(TRUE, nrow(result_table))
  }

  status <- rep(if (has_padj) "Not significant" else "Below cutoff", nrow(result_table))
  status[significant & !is.na(result_table$log2FoldChange) & result_table$log2FoldChange > lfc_cutoff] <- "Up regulated"
  status[significant & !is.na(result_table$log2FoldChange) & result_table$log2FoldChange < -lfc_cutoff] <- "Down regulated"

  levels <- c("Down regulated", if (has_padj) "Not significant" else "Below cutoff", "Up regulated")
  counts <- table(factor(status, levels = levels))
  data.frame(
    direction = names(counts),
    genes = as.integer(counts),
    stringsAsFactors = FALSE
  )
}

make_regulation_summary_plot <- function(result_table, alpha = 0.05, lfc_cutoff = 1) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }

  summary_data <- regulation_summary_counts(result_table, alpha, lfc_cutoff)
  has_padj <- "padj" %in% names(result_table) && any(!is.na(result_table$padj))
  fold_change_cutoff <- signif(2^lfc_cutoff, 4)
  subtitle <- if (has_padj) {
    sprintf("Criteria: adjusted p < %s and fold-change > %s", alpha, fold_change_cutoff)
  } else {
    sprintf("Criteria: fold-change > %s; no adjusted p-values supplied", fold_change_cutoff)
  }

  ggplot2::ggplot(
    summary_data,
    ggplot2::aes(x = direction, y = genes, fill = direction)
  ) +
    ggplot2::geom_col(width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = genes), vjust = -0.35, size = 4) +
    ggplot2::scale_fill_manual(
      values = c(
        "Down regulated" = "#2563EB",
        "Not significant" = "#6B7280",
        "Below cutoff" = "#6B7280",
        "Up regulated" = "#B42318"
      )
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      x = NULL,
      y = "Gene count",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11)
    )
}

format_result_preview <- function(result_table, n = 100) {
  preview <- utils::head(result_table, n)

  if ("description" %in% names(preview)) {
    preview$description <- ifelse(
      is.na(preview$description) | preview$description == "",
      "",
      ifelse(nchar(preview$description) > 120, paste0(substr(preview$description, 1, 117), "..."), preview$description)
    )
  }

  empty_model_cols <- intersect(c("lfcSE", "stat"), names(preview))
  empty_model_cols <- empty_model_cols[vapply(preview[empty_model_cols], function(column) all(is.na(column)), logical(1))]
  if (length(empty_model_cols) > 0) {
    preview <- preview[setdiff(names(preview), empty_model_cols)]
  }

  pvalue_cols <- intersect(c("pvalue", "padj"), names(preview))
  for (column_name in pvalue_cols) {
    preview[[column_name]] <- ifelse(
      is.na(preview[[column_name]]),
      "NA",
      formatC(preview[[column_name]], format = "e", digits = 3)
    )
  }

  numeric_cols <- vapply(preview, is.numeric, logical(1))
  preview[numeric_cols] <- lapply(preview[numeric_cols], function(column) signif(column, 4))
  preview
}

make_ma_like_plot <- function(workflow_result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is not installed. Run scripts/install_packages.R first.", call. = FALSE)
  }

  plot_data <- workflow_result$results
  plot_data$display_base_mean <- plot_data$baseMean
  if (!any(!is.na(plot_data$display_base_mean))) {
    plot_data$display_base_mean <- seq_len(nrow(plot_data))
  }

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = display_base_mean, y = log2FoldChange)
  ) +
    ggplot2::geom_point(alpha = 0.72, size = 1.8, color = "#0F766E", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#6B7280") +
    ggplot2::scale_x_continuous(trans = if (any(plot_data$display_base_mean > 0, na.rm = TRUE)) "log10" else "identity") +
    ggplot2::labs(
      x = if (any(!is.na(workflow_result$results$baseMean))) "Mean expression" else "Gene rank",
      y = "log2 fold change"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

draw_ma_plot <- function(workflow_result) {
  if (!is.null(workflow_result$deseq_result)) {
    DESeq2::plotMA(workflow_result$deseq_result, ylim = c(-5, 5), main = "MA plot")
  } else {
    print(make_ma_like_plot(workflow_result))
  }
}

write_result_bundle <- function(workflow_result, output_dir, lfc_cutoff = 1) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  utils::write.csv(
    workflow_result$results,
    file = file.path(output_dir, "deseq2_results.csv"),
    row.names = FALSE
  )
  if (!is.null(workflow_result$normalized_counts)) {
    utils::write.csv(
      workflow_result$normalized_counts,
      file = file.path(output_dir, "matrix_values.csv"),
      row.names = FALSE
    )
  }

  if (!is.null(workflow_result$pca_data)) {
    ggplot2::ggsave(
      filename = file.path(output_dir, "pca_plot.png"),
      plot = make_pca_plot(workflow_result),
      width = 7,
      height = 5,
      dpi = 160
    )
  }
  ggplot2::ggsave(
    filename = file.path(output_dir, "volcano_plot.png"),
    plot = make_volcano_plot(workflow_result$results, workflow_result$alpha, lfc_cutoff),
    width = 7,
    height = 5,
    dpi = 160
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "regulation_summary_plot.png"),
    plot = make_regulation_summary_plot(workflow_result$results, workflow_result$alpha, lfc_cutoff),
    width = 7,
    height = 5,
    dpi = 160
  )

  if (!is.null(workflow_result$deseq_result)) {
    grDevices::png(file.path(output_dir, "ma_plot.png"), width = 1120, height = 800, res = 160)
    on.exit(grDevices::dev.off(), add = TRUE)
    DESeq2::plotMA(workflow_result$deseq_result, ylim = c(-5, 5), main = "MA plot")
  } else {
    ggplot2::ggsave(
      filename = file.path(output_dir, "ma_plot.png"),
      plot = make_ma_like_plot(workflow_result),
      width = 7,
      height = 5,
      dpi = 160
    )
  }

  summary_lines <- c(
    "TranscriptoScope analysis summary",
    sprintf("Workflow: %s", if (is.null(workflow_result$workflow_label)) "DESeq2" else workflow_result$workflow_label),
    if (!is.null(workflow_result$reference_level)) sprintf("Reference group: %s", workflow_result$reference_level),
    if (!is.null(workflow_result$treatment_level)) sprintf("Comparison group: %s", workflow_result$treatment_level),
    sprintf("Genes tested: %s", workflow_result$genes_tested),
    if (!is.na(workflow_result$samples_tested)) sprintf("Samples tested: %s", workflow_result$samples_tested),
    sprintf("Adjusted p-value threshold: %s", workflow_result$alpha),
    if (length(workflow_result$warnings) > 0) {
      c("", "Warnings:", paste("- ", workflow_result$warnings))
    } else {
      c("", "Warnings: none")
    }
  )
  writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

  invisible(output_dir)
}
