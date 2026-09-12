#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else
  normalizePath(file.path(dirname(commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]), "..", "..", ".."), mustWork = FALSE)
out <- if (length(args) >= 2) args[[2]] else file.path(repo, "reports", "jmva_extension", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ordering_file <- file.path(repo, "reports", "jmva_extension", "ordering", "final", "ordering_summary.csv")
pe_file <- file.path(repo, "reports", "jmva_extension", "permutation_entropy", "final", "permutation_entropy_summary.csv")

ord <- read.csv(ordering_file, check.names = FALSE)
case_labels <- c(
  gaussian_rectangle_d5 = "Gaussian rectangle",
  gaussian_strong_d5 = "Strong Gaussian CDF",
  genuine_rvine_d5 = "Genuine R-vine"
)

pdf(file.path(out, "fig6_ordering_experiment.pdf"), width = 7.4, height = 4.4, useDingbats = FALSE)
par(mfrow = c(1, 3), mar = c(4.4, 4.2, 2.3, 0.7), oma = c(0, 0, 0, 0))
for (case in names(case_labels)) {
  d <- ord[ord$case_id == case, ]
  rel_v <- d$V_pi / min(d$V_pi)
  rel_r <- d$rqmc_rmse / min(d$rqmc_rmse)
  ymax <- max(c(rel_v, rel_r)) * 1.08
  x <- seq_len(nrow(d))
  plot(x, rel_v, type = "b", pch = 16, lty = 1, ylim = c(0.9, ymax),
       xaxt = "n", xlab = "Admissible order", ylab = "Ratio to best",
       main = case_labels[[case]], cex.main = 0.92)
  lines(x, rel_r, type = "b", pch = 1, lty = 2)
  axis(1, at = x, labels = d$order, las = 2, cex.axis = 0.78)
  abline(h = 1, col = "grey70", lty = 3)
  if (case == names(case_labels)[1])
    legend("topleft", c(expression(V[pi]), "RQMC RMSE"), pch = c(16, 1),
           lty = c(1, 2), bty = "n", cex = 0.82)
}
dev.off()

pe <- read.csv(pe_file, check.names = FALSE)
industries <- unique(pe$industry)
cols <- ifelse(pe$transform == "returns", "#2266AA", "#CC6677")
pch <- ifelse(pe$transform == "returns", 16, 17)

pdf(file.path(out, "fig7_permutation_entropy.pdf"), width = 7.4, height = 4.3, useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(4.2, 4.3, 2.1, 0.8))
lims <- range(c(pe$permutation_entropy_empirical, pe$permutation_entropy_model))
plot(pe$permutation_entropy_empirical, pe$permutation_entropy_model,
     col = cols, pch = pch, xlim = lims, ylim = lims,
     xlab = "Empirical normalized entropy", ylab = "Vine-model normalized entropy",
     main = "Permutation entropy")
abline(0, 1, col = "grey45", lty = 2)
legend("bottomright", c("Returns", "Absolute returns"), pch = c(16, 17),
       col = c("#2266AA", "#CC6677"), lty = c(1, 2), bty = "n",
       cex = 0.80, y.intersp = 0.8)

x <- seq_along(industries)
ret <- pe$total_variation[match(paste(industries, "returns"), paste(pe$industry, pe$transform))]
abv <- pe$total_variation[match(paste(industries, "absolute"), paste(pe$industry, pe$transform))]
plot(x, ret, type = "b", pch = 16, col = "#2266AA", xaxt = "n",
     ylim = range(c(ret, abv)), xlab = "Industry", ylab = "Total-variation distance",
     main = "Pattern-distribution discrepancy")
lines(x, abv, type = "b", pch = 17, col = "#CC6677", lty = 2)
axis(1, at = x, labels = industries, las = 2, cex.axis = 0.75)
dev.off()

cat("JMVA figures written to", normalizePath(out), "\n")
