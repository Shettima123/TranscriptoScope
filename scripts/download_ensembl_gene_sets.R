args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

gene_set_dir <- file.path(app_dir, "gene_sets")
dir.create(gene_set_dir, showWarnings = FALSE, recursive = TRUE)

species <- list(
  list(
    key = "yeast_scerevisiae_go",
    annotation_key = "yeast_scerevisiae",
    label = "Yeast GO terms (Saccharomyces cerevisiae)",
    scientific_name = "Saccharomyces cerevisiae",
    dataset = "scerevisiae_gene_ensembl",
    taxonomy_id = "559292",
    assembly = "R64-1-1",
    file = "yeast_scerevisiae_go.csv"
  ),
  list(
    key = "human_hsapiens_go",
    annotation_key = "human_hsapiens",
    label = "Human GO terms (Homo sapiens)",
    scientific_name = "Homo sapiens",
    dataset = "hsapiens_gene_ensembl",
    taxonomy_id = "9606",
    assembly = "GRCh38.p14",
    file = "human_hsapiens_go.csv"
  ),
  list(
    key = "fruitfly_dmelanogaster_go",
    annotation_key = "fruitfly_dmelanogaster",
    label = "Fruit fly GO terms (Drosophila melanogaster)",
    scientific_name = "Drosophila melanogaster",
    dataset = "dmelanogaster_gene_ensembl",
    taxonomy_id = "7227",
    assembly = "BDGP6.54",
    file = "fruitfly_dmelanogaster_go.csv"
  )
)

source_url <- "https://www.ensembl.org/biomart/martservice"
attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "go_id",
  "name_1006",
  "namespace_1003"
)

build_query <- function(dataset) {
  attribute_xml <- paste(sprintf('    <Attribute name="%s" />', attributes), collapse = "\n")
  sprintf(
    '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="%s" interface="default">
%s
  </Dataset>
</Query>',
    dataset,
    attribute_xml
  )
}

download_gene_sets <- function(entry) {
  query_file <- tempfile(fileext = ".xml")
  output_tsv <- tempfile(fileext = ".tsv")
  writeLines(build_query(entry$dataset), query_file, useBytes = TRUE)

  curl <- Sys.which("curl")
  if (!nzchar(curl)) {
    stop("curl was not found on PATH. Windows 10+ normally includes curl.exe.", call. = FALSE)
  }

  status <- system2(
    curl,
    args = c(
      "-L",
      "--fail",
      "--silent",
      "--show-error",
      "-X",
      "POST",
      "--data-urlencode",
      paste0("query@", query_file),
      source_url,
      "-o",
      output_tsv
    )
  )
  if (!identical(status, 0L)) {
    stop(sprintf("BioMart GO download failed for %s.", entry$dataset), call. = FALSE)
  }

  first_line <- readLines(output_tsv, n = 1, warn = FALSE)
  if (length(first_line) == 0 || grepl("^Query ERROR", first_line)) {
    stop(sprintf("BioMart returned an error for %s: %s", entry$dataset, first_line), call. = FALSE)
  }

  gene_sets <- utils::read.delim(
    output_tsv,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(gene_sets) <- c("gene_id", "gene_symbol", "term_id", "term_name", "domain")
  gene_sets$gene_id <- trimws(gene_sets$gene_id)
  gene_sets$gene_symbol <- trimws(gene_sets$gene_symbol)
  gene_sets$term_id <- trimws(gene_sets$term_id)
  gene_sets$term_name <- trimws(gene_sets$term_name)
  gene_sets$domain <- trimws(gene_sets$domain)
  gene_sets <- gene_sets[
    !is.na(gene_sets$gene_id) & gene_sets$gene_id != "" &
      !is.na(gene_sets$term_id) & gene_sets$term_id != "",
    ,
    drop = FALSE
  ]
  gene_sets$term_name[is.na(gene_sets$term_name) | gene_sets$term_name == ""] <- gene_sets$term_id[is.na(gene_sets$term_name) | gene_sets$term_name == ""]
  gene_sets <- unique(gene_sets)
  gene_sets <- gene_sets[order(gene_sets$domain, gene_sets$term_id, gene_sets$gene_id), , drop = FALSE]

  output_file <- file.path(gene_set_dir, entry$file)
  utils::write.csv(gene_sets, output_file, row.names = FALSE, na = "")

  tested_terms <- length(unique(gene_sets$term_id))
  message(sprintf("Wrote %s rows across %s terms to %s", nrow(gene_sets), tested_terms, output_file))
  c(rows = nrow(gene_sets), terms = tested_terms)
}

counts <- lapply(species, download_gene_sets)

manifest <- do.call(rbind, lapply(seq_along(species), function(index) {
  entry <- species[[index]]
  data.frame(
    key = entry$key,
    annotation_key = entry$annotation_key,
    label = entry$label,
    scientific_name = entry$scientific_name,
    dataset = entry$dataset,
    taxonomy_id = entry$taxonomy_id,
    assembly = entry$assembly,
    file = entry$file,
    source = source_url,
    downloaded_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC"),
    rows = counts[[index]][["rows"]],
    terms = counts[[index]][["terms"]],
    collection = "go",
    kegg_organism = "",
    stringsAsFactors = FALSE
  )
}))

existing_manifest_file <- file.path(gene_set_dir, "manifest.csv")
if (file.exists(existing_manifest_file)) {
  existing_manifest <- utils::read.csv(existing_manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  if ("collection" %in% names(existing_manifest)) {
    non_go_rows <- existing_manifest[existing_manifest$collection != "go", , drop = FALSE]
    if (nrow(non_go_rows) > 0) {
      manifest <- rbind(manifest, non_go_rows[names(manifest)])
    }
  }
}

utils::write.csv(manifest, file.path(gene_set_dir, "manifest.csv"), row.names = FALSE, na = "")
message(sprintf("Wrote manifest to %s", file.path(gene_set_dir, "manifest.csv")))
