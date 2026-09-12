suppressPackageStartupMessages(library(rvinecopulib))

out <- Sys.getenv("OUT", unset = "")
if (out == "") stop("Set OUT to the publication benchmark run directory")
manifest_path <- file.path(out, "case_manifest.csv")
references_path <- file.path(out, "references.csv")
if (!file.exists(manifest_path)) stop("Missing case_manifest.csv; run Julia publication setup first")
if (!file.exists(references_path)) stop("Missing references.csv; generate references first")

powers <- as.integer(strsplit(Sys.getenv("POWERS", unset = "8,10,12,14,16"), ",")[[1]])
reps <- as.integer(Sys.getenv("EXTERNAL_REPS", unset = "5"))
base_seed <- as.integer(Sys.getenv("EXTERNAL_BASE_SEED", unset = "90210"))

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
refs <- read.csv(references_path, stringsAsFactors = FALSE)

parse_vec <- function(x) as.numeric(strsplit(x, "\\|")[[1]])

make_gaussian_dvine <- function(d, rho) {
  pcs <- lapply(seq_len(d - 1), function(k) {
    rho_k <- rho / (1 + 0.12 * (k - 1))
    replicate(d - k, bicop_dist("gaussian", parameters = rho_k), simplify = FALSE)
  })
  vinecop_dist(pair_copulas = pcs, structure = dvine_structure(order = seq_len(d), trunc_lvl = d - 1))
}

mixed_pair <- function(index) {
  switch(
    as.character(index),
    "1" = bicop_dist("gaussian", parameters = 0.45),
    "2" = bicop_dist("clayton", parameters = 1.4),
    "3" = bicop_dist("gumbel", parameters = 1.3),
    "4" = bicop_dist("frank", parameters = 2.5),
    "5" = bicop_dist("joe", parameters = 1.35),
    stop("invalid mixed family index")
  )
}

make_mixed_dvine <- function(d) {
  pcs <- lapply(seq_len(d - 1), function(k) {
    level <- lapply(seq_len(d - k), function(i) {
      idx <- ((3 * k + i - 1) %% 5) + 1
      mixed_pair(idx)
    })
    if (k == d - 1) level <- list(bicop_dist("frank", parameters = 2.5))
    level
  })
  vinecop_dist(pair_copulas = pcs, structure = dvine_structure(order = seq_len(d), trunc_lvl = d - 1))
}

build_model <- function(row) {
  if (row$model != "D") return(NULL)
  if (row$design == "gaussian_moderate") return(make_gaussian_dvine(row$d, 0.45))
  if (row$design == "gaussian_strong") return(make_gaussian_dvine(row$d, 0.82))
  if (row$design == "mixed") return(make_mixed_dvine(row$d))
  NULL
}

rows <- list()
idx <- 1L
for (j in seq_len(nrow(manifest))) {
  case <- manifest[j, ]
  lower <- parse_vec(case$lower)
  upper <- parse_vec(case$upper)
  # pvinecop evaluates CDFs; general hyperrectangles are intentionally not
  # reconstructed by 2^d inclusion-exclusion.
  if (any(lower != 0)) next
  vc <- build_model(case)
  if (is.null(vc)) next
  ref <- refs[refs$case_id == case$case_id, ]
  if (nrow(ref) != 1) next

  # JIT/library warm-up outside recorded timing.
  invisible(pvinecop(matrix(upper, nrow = 1), vc, n_mc = 64, cores = 1))

  for (power in powers) {
    N <- 2^power
    for (rep in seq_len(reps)) {
      seed <- base_seed + rep - 1L
      set.seed(seed)
      gc()
      t0 <- proc.time()[["elapsed"]]
      estimate <- pvinecop(matrix(upper, nrow = 1), vc, n_mc = N, cores = 1)
      elapsed <- proc.time()[["elapsed"]] - t0
      rows[[idx]] <- data.frame(
        case_id = case$case_id,
        d = case$d,
        power = power,
        N = N,
        rep = rep,
        seed = seed,
        method = "rvinecopulib_pvinecop",
        estimate = as.numeric(estimate),
        reference = ref$estimate,
        abs_error = abs(as.numeric(estimate) - ref$estimate),
        runtime_s = elapsed,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  cat(sprintf("rvinecopulib: %s complete\n", case$case_id))
}

if (length(rows) == 0) stop("No supported CDF cases found")
result <- do.call(rbind, rows)
write.csv(result, file.path(out, "external_rvinecopulib.csv"), row.names = FALSE)

metadata <- c(
  paste("R:", R.version.string),
  paste("rvinecopulib:", as.character(packageVersion("rvinecopulib"))),
  paste("cores: 1"),
  paste("n_mc budgets:", paste(2^powers, collapse = ","))
)
writeLines(metadata, file.path(out, "external_rvinecopulib_environment.txt"))
cat("External rvinecopulib benchmark written under", out, "\n")
