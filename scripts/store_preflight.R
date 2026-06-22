args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)

required_packages <- c(
  "shiny",
  "ggplot2",
  "DESeq2",
  "SummarizedExperiment",
  "S4Vectors",
  "fgsea"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop(
    sprintf("Missing bundled package(s): %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

setwd(app_dir)
source(file.path(app_dir, "R", "deseq_helpers.R"), local = TRUE)

versions <- vapply(
  required_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
)

cat("TranscriptoScope Store preflight OK\n")
cat(sprintf("App directory: %s\n", app_dir))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("R library paths: %s\n", paste(.libPaths(), collapse = "; ")))
cat("Bundled package versions:\n")
for (pkg in names(versions)) {
  cat(sprintf("- %s %s\n", pkg, versions[[pkg]]))
}
