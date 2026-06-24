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

msigdb_release <- "2026.1.Hs"
gtrd_url <- sprintf(
  "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/%s/c3.tft.gtrd.v%s.symbols.gmt",
  msigdb_release,
  msigdb_release
)
output_file <- file.path(gene_set_dir, "human_hsapiens_gtrd.csv")
gmt_file <- tempfile(fileext = ".gmt")

download_gmt <- function(url, destfile) {
  curl <- Sys.which("curl")
  if (nzchar(curl)) {
    curl_args <- c(
      "--ssl-no-revoke",
      "-L",
      "--fail",
      "--silent",
      "--show-error",
      "--connect-timeout",
      "20",
      "--max-time",
      "120",
      "--retry",
      "3",
      "--retry-delay",
      "2",
      url,
      "-o",
      destfile
    )
    status <- system2(curl, args = curl_args)
    if (identical(status, 0L) && file.exists(destfile) && file.info(destfile)$size > 0) {
      return(invisible(TRUE))
    }
    warning("curl download failed; retrying with R download.file().", call. = FALSE)
  }

  utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
  if (!file.exists(destfile) || file.info(destfile)$size == 0) {
    stop("MSigDB GTRD download did not create a nonempty GMT file.", call. = FALSE)
  }
  invisible(TRUE)
}

read_gtrd_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  rows <- vector("list", length(lines))
  row_index <- 0L
  for (line in lines) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) {
      next
    }
    term_id <- trimws(fields[1])
    genes <- unique(trimws(fields[-c(1, 2)]))
    genes <- genes[!is.na(genes) & genes != ""]
    if (!nzchar(term_id) || length(genes) == 0) {
      next
    }
    row_index <- row_index + 1L
    rows[[row_index]] <- data.frame(
      gene_id = genes,
      gene_symbol = genes,
      term_id = term_id,
      term_name = term_id,
      domain = "tf_target_gtrd",
      stringsAsFactors = FALSE
    )
  }

  rows <- rows[seq_len(row_index)]
  if (length(rows) == 0) {
    stop("No valid TF.Target.GTRD gene sets were found in the GMT file.", call. = FALSE)
  }
  gene_sets <- unique(do.call(rbind, rows))
  gene_sets <- gene_sets[order(gene_sets$term_id, gene_sets$gene_id), , drop = FALSE]
  rownames(gene_sets) <- NULL
  gene_sets
}

download_gmt(gtrd_url, gmt_file)
gene_sets <- read_gtrd_gmt(gmt_file)
utils::write.csv(gene_sets, output_file, row.names = FALSE, na = "")

manifest_file <- file.path(gene_set_dir, "manifest.csv")
manifest <- if (file.exists(manifest_file)) {
  utils::read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}

gtrd_row <- data.frame(
  key = "human_hsapiens_gtrd",
  annotation_key = "human_hsapiens",
  label = "Human TF.Target.GTRD transcription factor targets (Homo sapiens)",
  scientific_name = "Homo sapiens",
  dataset = sprintf("MSigDB %s C3:TFT GTRD", msigdb_release),
  taxonomy_id = "9606",
  assembly = "GRCh38.p14",
  file = "human_hsapiens_gtrd.csv",
  source = gtrd_url,
  downloaded_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%d %H:%M:%S UTC"),
  rows = nrow(gene_sets),
  terms = length(unique(gene_sets$term_id)),
  collection = "gtrd",
  kegg_organism = "",
  stringsAsFactors = FALSE
)

if (nrow(manifest) > 0) {
  missing_cols <- setdiff(names(gtrd_row), names(manifest))
  for (col in missing_cols) {
    manifest[[col]] <- ""
  }
  extra_cols <- setdiff(names(manifest), names(gtrd_row))
  for (col in extra_cols) {
    gtrd_row[[col]] <- ""
  }
  manifest <- manifest[manifest$key != gtrd_row$key, names(manifest), drop = FALSE]
  manifest <- rbind(manifest, gtrd_row[names(manifest)])
} else {
  manifest <- gtrd_row
}

utils::write.csv(manifest, manifest_file, row.names = FALSE, na = "")
message(sprintf(
  "Wrote %s rows across %s TF.Target.GTRD terms to %s",
  nrow(gene_sets),
  length(unique(gene_sets$term_id)),
  output_file
))
message(sprintf("Updated manifest at %s", manifest_file))
