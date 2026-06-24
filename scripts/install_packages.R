setup_user_library <- function() {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) {
    version_dir <- paste(R.version$major, sub("\\..*$", "", R.version$minor), sep = ".")
    user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", version_dir)
  }

  user_lib <- path.expand(user_lib)
  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }

  if (!dir.exists(user_lib)) {
    stop(sprintf("Could not create per-user R library: %s", user_lib), call. = FALSE)
  }

  .libPaths(unique(c(normalizePath(user_lib, winslash = "/", mustWork = TRUE), .libPaths())))
  message(sprintf("Using per-user R library: %s", .libPaths()[1]))
}

setup_user_library()

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  install.packages.check.source = "no",
  install.packages.compile.from.source = "never"
)

ncpus <- max(1L, parallel::detectCores(logical = FALSE) - 1L)

cran_packages <- c("shiny", "ggplot2", "WGCNA")
bioc_packages <- c("DESeq2", "SummarizedExperiment", "S4Vectors", "fgsea", "impute", "preprocessCore")

install_missing_cran <- function(packages) {
  installed <- rownames(utils::installed.packages())
  missing <- setdiff(packages, installed)
  if (length(missing) > 0) {
    message(sprintf("Installing CRAN packages: %s", paste(missing, collapse = ", ")))
    utils::install.packages(missing, type = "binary", Ncpus = ncpus)
  } else {
    message("CRAN packages are already installed.")
  }
}

install_missing_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager.")
    utils::install.packages("BiocManager", type = "binary", Ncpus = ncpus)
  }
  installed <- rownames(utils::installed.packages())
  missing <- setdiff(packages, installed)
  if (length(missing) > 0) {
    message(sprintf("Installing Bioconductor packages: %s", paste(missing, collapse = ", ")))
    BiocManager::install(missing, ask = FALSE, update = FALSE, type = "binary", Ncpus = ncpus)
  } else {
    message("Bioconductor packages are already installed.")
  }
}

install_missing_bioc(bioc_packages)
install_missing_cran(cran_packages)

message("All required packages are installed.")
