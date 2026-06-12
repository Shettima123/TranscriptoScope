args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args) >= 1) as.integer(args[1]) else 7865

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

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
  app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  app_dir <- normalizePath(".", mustWork = TRUE)
}

local_lib <- normalizePath(file.path(app_dir, "..", "r-lib"), mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

shiny::runApp(app_dir, host = "127.0.0.1", port = port, launch.browser = FALSE)
