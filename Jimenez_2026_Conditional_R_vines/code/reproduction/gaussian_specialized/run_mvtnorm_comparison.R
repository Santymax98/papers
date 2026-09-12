if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("R package 'mvtnorm' is required")
}

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) != 1L) stop("cannot locate script directory")
here <- dirname(normalizePath(sub("^--file=", "", script_arg)))
out_dir <- Sys.getenv("OUT_DIR", unset = file.path(here, "..", "..", "..", "output", "gaussian_specialized"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

N <- 4096L
REPS <- 30L
BASE_SEED <- 829828L
ABSEPS <- 1e-6
RELEPS <- 1e-6

parse_vector <- function(x) as.numeric(strsplit(x, "|", fixed = TRUE)[[1L]])
inputs <- read.csv(file.path(out_dir, "mvtnorm_inputs.csv"), stringsAsFactors = FALSE)
raw <- vector("list", nrow(inputs) * REPS)
out_index <- 1L

for (case_index in seq_len(nrow(inputs))) {
  case <- inputs[case_index, ]
  d <- as.integer(case$d)
  lower <- parse_vector(case$lower_z)
  upper <- parse_vector(case$upper_z)
  sigma <- matrix(parse_vector(case$correlation), nrow = d, ncol = d)
  algorithm <- mvtnorm::GenzBretz(maxpts = N, abseps = ABSEPS, releps = RELEPS)

  set.seed(BASE_SEED)
  invisible(mvtnorm::pmvnorm(
    lower = lower, upper = upper, mean = numeric(d), sigma = sigma,
    algorithm = algorithm
  ))

  for (rep in seq_len(REPS)) {
    seed <- BASE_SEED + rep - 1L
    set.seed(seed)
    elapsed <- system.time({
      result <- mvtnorm::pmvnorm(
        lower = lower, upper = upper, mean = numeric(d), sigma = sigma,
        algorithm = algorithm
      )
    })[["elapsed"]]
    estimate <- as.numeric(result)
    internal_error <- as.numeric(attr(result, "error"))
    message <- as.character(attr(result, "msg"))
    reference <- as.numeric(case$reference)
    raw[[out_index]] <- data.frame(
      case_id = case$case_id,
      d = d,
      event = case$event,
      dependence = case$dependence,
      method = "mvtnorm_genzbretz",
      N_or_maxpts = N,
      replicate = rep,
      seed = seed,
      estimate = estimate,
      reference = reference,
      abs_error = abs(estimate - reference),
      runtime_s = elapsed,
      internal_error = internal_error,
      status = message,
      algorithm = "mvtnorm::GenzBretz",
      stringsAsFactors = FALSE
    )
    out_index <- out_index + 1L
  }
}

raw <- do.call(rbind, raw)
write.csv(raw, file.path(out_dir, "mvtnorm_raw_results.csv"), row.names = FALSE)

case_ids <- unique(raw$case_id)
summary_rows <- lapply(case_ids, function(case_id) {
  rows <- raw[raw$case_id == case_id, ]
  data.frame(
    case_id = case_id,
    d = rows$d[1L],
    event = rows$event[1L],
    dependence = rows$dependence[1L],
    method = rows$method[1L],
    N_or_maxpts = N,
    replicates = nrow(rows),
    reference = rows$reference[1L],
    rmse = sqrt(mean((rows$estimate - rows$reference)^2)),
    median_abs_error = median(rows$abs_error),
    median_runtime_s = median(rows$runtime_s),
    median_internal_error = median(rows$internal_error),
    tolerance = ABSEPS,
    successful_statuses = sum(rows$status == "Normal Completion"),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(out_dir, "mvtnorm_summary.csv"), row.names = FALSE)

cat("mvtnorm comparison complete\n")
cat("  R             =", as.character(getRversion()), "\n")
cat("  mvtnorm       =", as.character(utils::packageVersion("mvtnorm")), "\n")
cat("  cases         =", length(case_ids), "\n")
cat("  rows          =", nrow(raw), "\n")
cat("  normal status =", sum(raw$status == "Normal Completion"), "/", nrow(raw), "\n")

