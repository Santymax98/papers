required <- c(mvtnorm = "1.3-3", rvinecopulib = "0.7.3.1.0")

for (package in names(required)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("missing R package %s; run R/install_exact_versions.R", package))
  }
  installed <- packageDescription(package)$Version
  if (!identical(installed, required[[package]])) {
    stop(sprintf("%s version %s installed; required %s",
                 package, installed, required[[package]]))
  }
  cat(sprintf("%-14s %s\n", package, installed))
}
cat(sprintf("R              %s\n", as.character(getRversion())))
if (getRversion() != "4.5.3") {
  stop("R 4.5.3 is required for the recorded environment")
}
cat("R environment check passed.\n")

