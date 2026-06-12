args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

annotation_dir <- file.path(app_dir, "annotations")
dir.create(annotation_dir, showWarnings = FALSE, recursive = TRUE)

species <- list(
  list(
    key = "yeast_scerevisiae",
    label = "Yeast (Saccharomyces cerevisiae)",
    scientific_name = "Saccharomyces cerevisiae",
    dataset = "scerevisiae_gene_ensembl",
    taxonomy_id = "559292",
    assembly = "R64-1-1",
    file = "yeast_scerevisiae_ensembl.csv"
  ),
  list(
    key = "human_hsapiens",
    label = "Human (Homo sapiens)",
    scientific_name = "Homo sapiens",
    dataset = "hsapiens_gene_ensembl",
    taxonomy_id = "9606",
    assembly = "GRCh38.p14",
    file = "human_hsapiens_ensembl.csv"
  ),
  list(
    key = "fruitfly_dmelanogaster",
    label = "Fruit fly (Drosophila melanogaster)",
    scientific_name = "Drosophila melanogaster",
    dataset = "dmelanogaster_gene_ensembl",
    taxonomy_id = "7227",
    assembly = "BDGP6.54",
    file = "fruitfly_dmelanogaster_ensembl.csv"
  )
)

source_url <- "https://www.ensembl.org/biomart/martservice"
attributes <- c(
  "ensembl_gene_id",
  "external_gene_name",
  "description",
  "gene_biotype",
  "chromosome_name",
  "start_position",
  "end_position"
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

download_annotation <- function(entry) {
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
    stop(sprintf("BioMart download failed for %s.", entry$dataset), call. = FALSE)
  }

  first_line <- readLines(output_tsv, n = 1, warn = FALSE)
  if (length(first_line) == 0 || grepl("^Query ERROR", first_line)) {
    stop(sprintf("BioMart returned an error for %s: %s", entry$dataset, first_line), call. = FALSE)
  }

  annotation <- utils::read.delim(
    output_tsv,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(annotation) <- c(
    "ensembl_gene_id",
    "gene_symbol",
    "description",
    "gene_biotype",
    "chromosome",
    "start",
    "end"
  )
  annotation$ensembl_gene_id <- trimws(annotation$ensembl_gene_id)
  annotation$gene_symbol <- trimws(annotation$gene_symbol)
  annotation <- annotation[!is.na(annotation$ensembl_gene_id) & annotation$ensembl_gene_id != "", , drop = FALSE]
  annotation <- unique(annotation)
  annotation <- annotation[order(annotation$ensembl_gene_id), , drop = FALSE]

  output_file <- file.path(annotation_dir, entry$file)
  utils::write.csv(annotation, output_file, row.names = FALSE, na = "")
  message(sprintf("Wrote %s rows to %s", nrow(annotation), output_file))
  nrow(annotation)
}

row_counts <- vapply(species, download_annotation, integer(1))

manifest <- do.call(rbind, lapply(seq_along(species), function(index) {
  entry <- species[[index]]
  data.frame(
    key = entry$key,
    label = entry$label,
    scientific_name = entry$scientific_name,
    dataset = entry$dataset,
    taxonomy_id = entry$taxonomy_id,
    assembly = entry$assembly,
    file = entry$file,
    source = source_url,
    downloaded_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC"),
    rows = row_counts[[index]],
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(manifest, file.path(annotation_dir, "manifest.csv"), row.names = FALSE, na = "")
message(sprintf("Wrote manifest to %s", file.path(annotation_dir, "manifest.csv")))
