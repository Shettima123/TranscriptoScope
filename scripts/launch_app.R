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

required_packages <- c("shiny", "ggplot2", "fgsea")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required package(s): %s. Run scripts/install_packages.R first.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

setwd(app_dir)
shiny::runApp(app_dir, launch.browser = TRUE)
