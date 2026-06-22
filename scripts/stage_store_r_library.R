args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: stage_store_r_library.R <target-lib> <comma-separated-root-packages>", call. = FALSE)
}

target_lib <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
root_packages <- trimws(strsplit(args[[2]], ",", fixed = TRUE)[[1]])
root_packages <- root_packages[nzchar(root_packages)]

if (length(root_packages) == 0) {
  stop("At least one root package is required.", call. = FALSE)
}

installed <- utils::installed.packages()
dependencies <- tools::package_dependencies(
  root_packages,
  db = installed,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)
packages <- sort(unique(c(root_packages, unlist(dependencies, use.names = FALSE))))

missing_packages <- setdiff(packages, rownames(installed))
if (length(missing_packages) > 0) {
  stop(
    sprintf("Missing installed package(s): %s", paste(missing_packages, collapse = ", ")),
    call. = FALSE
  )
}

base_or_recommended <- rownames(installed)[installed[, "Priority"] %in% c("base", "recommended")]
packages_to_copy <- setdiff(packages, base_or_recommended)

if (dir.exists(target_lib)) {
  unlink(target_lib, recursive = TRUE, force = TRUE)
}
dir.create(target_lib, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(target_lib)) {
  stop(sprintf("Could not create target library: %s", target_lib), call. = FALSE)
}

manifest <- data.frame(
  Package = packages_to_copy,
  Version = installed[packages_to_copy, "Version"],
  SourceLibPath = installed[packages_to_copy, "LibPath"],
  SourcePackageDir = file.path(installed[packages_to_copy, "LibPath"], packages_to_copy),
  stringsAsFactors = FALSE,
  row.names = NULL
)

for (source_dir in manifest$SourcePackageDir) {
  if (!dir.exists(source_dir)) {
    stop(sprintf("Package directory does not exist: %s", source_dir), call. = FALSE)
  }
}

utils::write.csv(
  manifest,
  file.path(target_lib, "TRANSCRIPTOSCOPE_STORE_R_PACKAGES.csv"),
  row.names = FALSE
)

cat(sprintf("Wrote manifest for %d R package(s) in %s\n", nrow(manifest), target_lib))
