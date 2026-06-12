args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

source(file.path(app_dir, "R", "deseq_helpers.R"))

gene_set_dir <- file.path(app_dir, "gene_sets")
manifest <- read_gene_set_manifest(gene_set_dir)
kegg_rows <- manifest[manifest$collection == "kegg", , drop = FALSE]

if (nrow(kegg_rows) == 0) {
  stop("No KEGG gene set databases are listed in gene_sets/manifest.csv.", call. = FALSE)
}

for (key in kegg_rows$key) {
  message(sprintf("Downloading %s", key))
  gene_sets <- read_builtin_gene_sets(
    gene_set_dir = gene_set_dir,
    gene_set_key = key,
    domain = "kegg_pathway"
  )
  message(sprintf(
    "Cached %s rows across %s pathways for %s",
    nrow(gene_sets$table),
    length(unique(gene_sets$table$term_id)),
    key
  ))
}
