# Publication preflight for the external rvinecopulib comparison.
# Fails if the R construction is not density-equivalent to the Julia model.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  RCPP_PARALLEL_NUM_THREADS = "1"
)

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_path))

if (!requireNamespace("rvinecopulib", quietly = TRUE)) {
  stop("Package 'rvinecopulib' is not installed. Run install.packages('rvinecopulib').")
}
suppressPackageStartupMessages(library(rvinecopulib))
source(file.path(script_dir, "rvinecopulib_common.R"))

out <- Sys.getenv("OUT", unset = "")
if (out == "") stop("Set OUT to the completed publication benchmark directory")
ext_dir <- file.path(out, "external_rvinecopulib")
cases_path <- file.path(ext_dir, "cases.csv")
audit_path <- file.path(ext_dir, "julia_density_audit.csv")
if (!file.exists(cases_path)) stop("Missing external cases.csv; run prepare_rvinecopulib_external.jl first")
if (!file.exists(audit_path)) stop("Missing julia_density_audit.csv; run prepare_rvinecopulib_external.jl first")

cases <- read.csv(cases_path, stringsAsFactors = FALSE)
audit <- read.csv(audit_path, stringsAsFactors = FALSE)

atol <- as.numeric(Sys.getenv("EXTERNAL_MODEL_AUDIT_ATOL", unset = "1e-7"))
rtol <- as.numeric(Sys.getenv("EXTERNAL_MODEL_AUDIT_RTOL", unset = "1e-7"))
rows <- list()
idx <- 1L

for (j in seq_len(nrow(cases))) {
  case <- cases[j, ]
  vc <- build_publication_model(case)
  aa <- audit[audit$case_id == case$case_id, ]
  if (nrow(aa) != 2L) stop(sprintf("expected two Julia audit points for %s", case$case_id))

  for (k in seq_len(nrow(aa))) {
    u <- parse_vec(aa$u[k])
    r_logpdf <- log(as.numeric(dvinecop(matrix(u, nrow = 1L), vc, cores = 1L)))
    j_logpdf <- as.numeric(aa$julia_logpdf[k])
    abs_diff <- abs(r_logpdf - j_logpdf)
    tolerance <- atol + rtol * abs(j_logpdf)
    passed <- is.finite(r_logpdf) && abs_diff <= tolerance
    rows[[idx]] <- data.frame(
      case_id = case$case_id,
      d = case$d,
      audit_point = aa$audit_point[k],
      julia_logpdf = j_logpdf,
      rvinecopulib_logpdf = r_logpdf,
      abs_difference = abs_diff,
      tolerance = tolerance,
      passed = passed,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
    if (!passed) {
      stop(sprintf(
        "model parity failed for %s audit point %d: Julia=%.17g R=%.17g diff=%.3e tol=%.3e",
        case$case_id, aa$audit_point[k], j_logpdf, r_logpdf, abs_diff, tolerance
      ))
    }
  }
  cat(sprintf("[%d/%d] %-36s model parity OK\n", j, nrow(cases), case$case_id))
}

model_audit <- do.call(rbind, rows)
write.csv(model_audit, file.path(ext_dir, "model_parity_audit.csv"), row.names = FALSE)

# Reproducibility audit for pvinecop's randomized QMC path. pvinecop forwards
# package-generated seeds to the C++ implementation; set.seed() should therefore
# make an individual call reproducible.
case <- cases[1, ]
vc <- build_publication_model(case)
u <- matrix(parse_vec(case$upper), nrow = 1L)
seed <- 190887L
set.seed(seed)
p1 <- as.numeric(pvinecop(u, vc, n_mc = 256L, cores = 1L))
set.seed(seed)
p2 <- as.numeric(pvinecop(u, vc, n_mc = 256L, cores = 1L))
set.seed(seed + 1L)
p3 <- as.numeric(pvinecop(u, vc, n_mc = 256L, cores = 1L))
same_seed_identical <- identical(p1, p2)
if (!same_seed_identical) stop("pvinecop reproducibility audit failed: identical R seeds gave different estimates")

audit_lines <- c(
  sprintf("R: %s", R.version.string),
  sprintf("rvinecopulib: %s", as.character(packageVersion("rvinecopulib"))),
  sprintf("CDF cases: %d", nrow(cases)),
  sprintf("density audit rows: %d", nrow(model_audit)),
  sprintf("density parity: all passed (atol=%g, rtol=%g)", atol, rtol),
  sprintf("same-seed pvinecop reproducibility: %s", same_seed_identical),
  sprintf("different-seed estimate changed: %s", !identical(p1, p3)),
  sprintf("seed audit p1: %.17g", p1),
  sprintf("seed audit p3: %.17g", p3)
)
writeLines(audit_lines, file.path(ext_dir, "PREFLIGHT.txt"))
cat("rvinecopulib external preflight passed for all CDF cases.\n")
