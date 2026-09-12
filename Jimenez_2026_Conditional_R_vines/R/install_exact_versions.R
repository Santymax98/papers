options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

required <- c(mvtnorm = "1.3-3", rvinecopulib = "0.7.3.1.0")
for (package in names(required)) {
  installed <- if (requireNamespace(package, quietly = TRUE)) {
    packageDescription(package)$Version
  } else {
    ""
  }
  if (!identical(installed, required[[package]])) {
    remotes::install_version(package, version = required[[package]], upgrade = "never")
  }
}

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) != 1L) stop("cannot locate installer directory")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
source(file.path(script_dir, "check_versions.R"))
