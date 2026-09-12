# Shared construction helpers for the publication comparison with rvinecopulib.
# The model definitions mirror publication_common.jl exactly.

parse_vec <- function(x) as.numeric(strsplit(x, "\\|", fixed = FALSE)[[1]])

make_gaussian_dvine <- function(d, rho) {
  pcs <- lapply(seq_len(d - 1), function(k) {
    rho_k <- rho / (1 + 0.12 * (k - 1))
    replicate(
      d - k,
      bicop_dist("gaussian", parameters = rho_k),
      simplify = FALSE
    )
  })
  vinecop_dist(
    pair_copulas = pcs,
    structure = dvine_structure(order = seq_len(d), trunc_lvl = d - 1)
  )
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
      # Julia: pool[mod1(3k + i, length(pool))]
      idx <- ((3 * k + i - 1) %% 5) + 1
      mixed_pair(idx)
    })
    if (k == d - 1) level <- list(bicop_dist("frank", parameters = 2.5))
    level
  })
  vinecop_dist(
    pair_copulas = pcs,
    structure = dvine_structure(order = seq_len(d), trunc_lvl = d - 1)
  )
}

make_genuine_rvine5 <- function() {
  ord <- c(1L, 3L, 2L, 4L, 5L)
  S <- list(
    c(2L, 2L, 4L, 5L),
    c(3L, 4L, 5L),
    c(4L, 5L),
    c(5L)
  )
  pcs <- list(
    list(
      bicop_dist("gaussian", parameters = 0.35),
      bicop_dist("clayton", parameters = 1.4),
      bicop_dist("frank", parameters = 2.0),
      bicop_dist("gumbel", parameters = 1.25)
    ),
    list(
      bicop_dist("frank", parameters = 1.7),
      bicop_dist("gaussian", parameters = -0.25),
      bicop_dist("clayton", parameters = 1.1)
    ),
    list(
      bicop_dist("gumbel", parameters = 1.3),
      bicop_dist("gaussian", parameters = 0.20)
    ),
    list(bicop_dist("clayton", parameters = 0.8))
  )
  vinecop_dist(
    pair_copulas = pcs,
    structure = rvine_structure(order = ord, struct_array = S)
  )
}

build_publication_model <- function(case) {
  if (case$model == "D" && case$design == "gaussian_moderate") {
    return(make_gaussian_dvine(case$d, 0.45))
  }
  if (case$model == "D" && case$design == "gaussian_strong") {
    return(make_gaussian_dvine(case$d, 0.82))
  }
  if (case$model == "D" && case$design == "mixed") {
    return(make_mixed_dvine(case$d))
  }
  if (case$model == "R" && case$design == "genuine_rvine" && case$d == 5) {
    return(make_genuine_rvine5())
  }
  stop(sprintf(
    "unsupported external case: id=%s model=%s design=%s d=%d",
    case$case_id, case$model, case$design, case$d
  ))
}

external_seed <- function(base_seed, case_index, power, rep) {
  # Unique, deterministic seed for every external randomized-QMC call.
  as.integer(base_seed + 100000L * case_index + 1000L * power + rep)
}

append_external_row <- function(path, row) {
  write.table(
    row,
    file = path,
    append = file.exists(path),
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    quote = TRUE,
    qmethod = "double"
  )
}

completed_external_keys <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- read.csv(path, stringsAsFactors = FALSE)
  if (nrow(x) == 0) return(character(0))
  paste(x$case_id, x$N, x$rep, sep = "|")
}
