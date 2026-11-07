# =============================================================================
# reproduce.R  (bundled with the multistrata app; not attached on launch)
#
# Optional utilities to (a) regenerate the two example CSVs shipped with the
# app and (b) reproduce the bivariate Khan & Reddy design (Case 1 / Case 2 and
# its data-free analogue) and print the package solution next to the paper's
# LINGO results.  Source this file from an R session to use the functions:
#
#     source(system.file("shiny", "reproduce.R", package = "multistrata"))
#     reproduce_khan_reddy()
# =============================================================================

# Regenerate the two example datasets used in the app's upload demo.
write_khan_reddy_data <- function(dir = ".", N = 2000L, seed = 20220502L) {
  set.seed(seed)
  x <- stats::rexp(N, rate = 1); x[x > 10] <- 10
  y <- 2 - 2 * sqrt(1 - stats::runif(N))
  z <- 0.5 * x + 1.5 * y + stats::rnorm(N, 0, 1)
  x <- round(x, 6); y <- round(y, 6); z <- round(z, 6)
  case1 <- data.frame(x = x, y = y, x_study = x, y_study = y)  # X={x,y}, Y={x_study,y_study}
  case2 <- data.frame(x = x, y = y, z = z)                     # X={x,y}, Y={z}
  f1 <- file.path(dir, "khan_reddy_case1.csv")
  f2 <- file.path(dir, "khan_reddy_case2.csv")
  utils::write.csv(case1, f1, row.names = FALSE)
  utils::write.csv(case2, f2, row.names = FALSE)
  message("Wrote example datasets:\n  ", f1, "\n  ", f2)
  invisible(c(case1 = f1, case2 = f2))
}

# Reproduce the paper's two data-based cases through the package.
reproduce_khan_reddy <- function(N = 2000L, n_total = 700L, n_starts = 30L,
                                 seed = 20220502L) {
  stopifnot(requireNamespace("multistrata", quietly = TRUE))
  set.seed(seed)
  x <- stats::rexp(N, rate = 1); x[x > 10] <- 10
  y <- 2 - 2 * sqrt(1 - stats::runif(N))
  z <- 0.5 * x + 1.5 * y + stats::rnorm(N, 0, 1)
  pop <- data.frame(x = x, y = y, z = z)
  dl  <- list(x = list(family = "exp",    rate = 1),
              y = list(family = "triang", a = 0, b = 2, c = 0))
  H_vec <- c(4L, 3L)

  unit_sum <- function(res) {
    ss <- res$stratum_stats
    s2 <- as.matrix(ss[, c("sigma2_k1", "sigma2_k2")])
    multistrata::composite_variance(s2, ss$W_h, res$allocation, ss$N_h,
                                       omega = c(1, 1))
  }

  cat("\n----- CASE 1: x and y are the study variables -----\n")
  c1 <- multistrata::stratify_multi(
    pop[, c("x", "y")], pop[, c("x", "y")], H_vec = H_vec, n_total = n_total,
    omega = c(0.5, 0.5), n_starts = n_starts, distributions = dl)
  cat(sprintf("  var(x_bar)+var(y_bar) = %.8f   [paper: 0.00003691]\n", unit_sum(c1)))
  cat("  x cuts:", paste(sprintf("%.4f", c1$boundaries[[1]]), collapse = ", "),
      "  [paper: 0.3917, 0.8938, 1.6092]\n")
  cat("  y cuts:", paste(sprintf("%.4f", c1$boundaries[[2]]), collapse = ", "),
      "  [paper: 0.3088, 0.7834]\n")

  cat("\n----- CASE 2: study variable z = 0.5x + 1.5y + e -----\n")
  c2 <- multistrata::stratify_multi(
    pop[, c("x", "y")], data.frame(z = z), H_vec = H_vec, n_total = n_total,
    omega = 1, n_starts = n_starts, distributions = dl)
  cat(sprintf("  var(z_bar) = %.8f   [paper: 0.00095293]\n", c2$objective))
  cat("  x cuts:", paste(sprintf("%.4f", c2$boundaries[[1]]), collapse = ", "),
      "  [paper: 0.7422, 1.7344, 3.2760]\n")
  cat("  y cuts:", paste(sprintf("%.4f", c2$boundaries[[2]]), collapse = ", "),
      "  [paper: 0.6021, 1.6227]\n\n")
  invisible(list(population = pop, case1 = c1, case2 = c2))
}

# Data-free analogue of Case 1 (density mode only, no data).
reproduce_khan_reddy_density <- function(n_total = 700L, N = 2000L,
                                         n_starts = 20L) {
  stopifnot(requireNamespace("multistrata", quietly = TRUE))
  res <- multistrata::stratify_multi_density(
    distributions = list(x = list(family = "exp",    rate = 1),
                         y = list(family = "triang", a = 0, b = 2, c = 0)),
    H_vec   = c(4L, 3L), n_total = n_total, N = N,
    ranges  = list(c(0, 10), c(0, 2)),
    omega   = c(0.5, 0.5), n_starts = n_starts, verbose = FALSE)
  cat("\n----- CASE 1, DATA-FREE (density mode) -----\n")
  cat("  x cuts:", paste(sprintf("%.4f", res$boundaries[[1]]), collapse = ", "),
      "  [paper: 0.3917, 0.8938, 1.6092]\n")
  cat("  y cuts:", paste(sprintf("%.4f", res$boundaries[[2]]), collapse = ", "),
      "  [paper: 0.3088, 0.7834]\n")
  cat(sprintf("  Composite Phi (omega .5/.5) = %.8f\n", res$objective))
  invisible(res)
}
