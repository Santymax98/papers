# Final external publication benchmark against rvinecopulib::pvinecop().
# CDF cases only.  Model construction and warm-up are excluded from timing.

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
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)

cases_path <- file.path(ext_dir, "cases.csv")
preflight_path <- file.path(ext_dir, "model_parity_audit.csv")
if (!file.exists(cases_path)) stop("Missing cases.csv; run prepare_rvinecopulib_external.jl first")
if (!file.exists(preflight_path)) stop("Missing model_parity_audit.csv; run rvinecopulib_preflight.R first")

preflight <- read.csv(preflight_path, stringsAsFactors = FALSE)
if (!all(preflight$passed)) stop("Model parity preflight contains failures")

powers <- as.integer(strsplit(Sys.getenv("POWERS", unset = "8,10,12,14,16"), ",")[[1]])
reps <- as.integer(Sys.getenv("EXTERNAL_REPS", unset = "30"))
base_seed <- as.integer(Sys.getenv("EXTERNAL_BASE_SEED", unset = "90210"))
resume <- tolower(Sys.getenv("EXTERNAL_RESUME", unset = "true")) %in% c("1", "true", "t", "yes", "y")
cores <- as.integer(Sys.getenv("EXTERNAL_CORES", unset = "1"))
if (cores != 1L) stop("Publication external timing must use EXTERNAL_CORES=1")

cases <- read.csv(cases_path, stringsAsFactors = FALSE)
raw_path <- file.path(ext_dir, "raw_results.csv")
if (!resume && file.exists(raw_path)) file.remove(raw_path)
done <- if (resume) completed_external_keys(raw_path) else character(0)

cat("External rvinecopulib publication benchmark\n")
cat("  cases      =", nrow(cases), "\n")
cat("  powers     =", paste(powers, collapse = ","), "\n")
cat("  reps       =", reps, "\n")
cat("  cores      =", cores, "\n")
cat("  output     =", ext_dir, "\n")
cat("  package    =", as.character(packageVersion("rvinecopulib")), "\n")

execution_index <- if (file.exists(raw_path)) nrow(read.csv(raw_path, stringsAsFactors = FALSE)) else 0L

for (j in seq_len(nrow(cases))) {
  case <- cases[j, ]
  vc <- build_publication_model(case)
  u <- matrix(parse_vec(case$upper), nrow = 1L)

  # Compiled-library/dynamic-loader/cache warm-up, deliberately not recorded.
  set.seed(external_seed(base_seed, j, 6L, 0L))
  invisible(pvinecop(u, vc, n_mc = 64L, cores = cores))

  cat(sprintf("[%d/%d] %s\n", j, nrow(cases), case$case_id))
  for (power in powers) {
    N <- as.integer(2^power)
    invisible(gc())
    for (rep in seq_len(reps)) {
      key <- paste(case$case_id, N, rep, sep = "|")
      if (key %in% done) next

      seed <- external_seed(base_seed, j, power, rep)
      set.seed(seed)
      t0 <- proc.time()[["elapsed"]]
      estimate <- as.numeric(pvinecop(u, vc, n_mc = N, cores = cores))
      elapsed <- proc.time()[["elapsed"]] - t0
      execution_index <- execution_index + 1L

      row <- data.frame(
        case_id = case$case_id,
        design = case$design,
        model = case$model,
        point = case$point,
        d = case$d,
        power = power,
        N = N,
        rep = rep,
        seed = seed,
        method = "rvinecopulib_pvinecop",
        execution_index = execution_index,
        estimate = estimate,
        reference = case$reference,
        reference_uncertainty = case$reference_uncertainty,
        reference_target_met = case$reference_target_met,
        abs_error = abs(estimate - case$reference),
        runtime_s = elapsed,
        cores = cores,
        stringsAsFactors = FALSE
      )
      append_external_row(raw_path, row)
      done <- c(done, key)
    }
    cat(sprintf("    N=%d complete\n", N))
    flush.console()
  }
}

# Environment and implementation provenance.  pvinecop's public documentation
# identifies n_mc as the number of quasi-Monte Carlo integration samples.
env_lines <- c(
  sprintf("created_at: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  sprintf("R: %s", R.version.string),
  sprintf("platform: %s", R.version$platform),
  sprintf("rvinecopulib: %s", as.character(packageVersion("rvinecopulib"))),
  sprintf("cores: %d", cores),
  sprintf("reps: %d", reps),
  sprintf("n_mc_budgets: %s", paste(2^powers, collapse = ",")),
  sprintf("base_seed: %d", base_seed),
  sprintf("RNGkind: %s", paste(RNGkind(), collapse = ",")),
  sprintf("OMP_NUM_THREADS: %s", Sys.getenv("OMP_NUM_THREADS")),
  sprintf("OPENBLAS_NUM_THREADS: %s", Sys.getenv("OPENBLAS_NUM_THREADS")),
  sprintf("VECLIB_MAXIMUM_THREADS: %s", Sys.getenv("VECLIB_MAXIMUM_THREADS")),
  sprintf("RCPP_PARALLEL_NUM_THREADS: %s", Sys.getenv("RCPP_PARALLEL_NUM_THREADS"))
)
writeLines(env_lines, file.path(ext_dir, "environment.txt"))
capture.output(sessionInfo(), file = file.path(ext_dir, "sessionInfo.txt"))
cat("External rvinecopulib raw benchmark complete.\n")
