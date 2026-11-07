# End-to-end tests for stratify_multi() -------------------------------------

make_pop <- function(N = 1500) {
  set.seed(31)
  sim_population(N = N, p = 2, q = 2, scenario = "skewed")
}

test_that("empirical mode returns a coherent solution", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  res <- stratify_multi(
    data = pop[, c("X1", "X2")], study_vars = pop[, c("Y1", "Y2")],
    H_vec = c(3, 3), n_total = 300, n_starts = 3, max_iter = 25)

  expect_s3_class(res, "multistrata")
  expect_equal(sum(res$allocation), 300)
  expect_length(res$boundaries, 2)
  for (j in 1:2) {
    b <- res$boundaries[[j]]
    expect_length(b, 2)
    expect_true(all(diff(b) > 0))
    expect_true(all(b > min(pop[, j]) & b < max(pop[, j])))
  }
  expect_true(is.finite(res$objective) && res$objective > 0)
  expect_true(is.finite(res$optimality_gap) && res$optimality_gap >= 0)
  expect_true(is.finite(res$kkt_residual))
  # joint optimum should improve on naive equal allocation over same strata
  ss   <- res$stratum_stats
  s2   <- as.matrix(ss[, c("sigma2_k1", "sigma2_k2")])
  neq  <- rep(300 / nrow(ss), nrow(ss))
  expect_lte(res$objective,
             composite_variance(s2, ss$W_h, neq, ss$N_h, c(0.5, 0.5)) + 1e-12)
})

test_that("separate KKT residuals and per-start log are reported", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  res <- stratify_multi(
    data = pop[, c("X1", "X2")], study_vars = pop[, c("Y1", "Y2")],
    H_vec = c(3, 3), n_total = 300, n_starts = 4, max_iter = 25)

  # allocation and boundary residuals are present and non-negative
  expect_true(is.finite(res$kkt_alloc) && res$kkt_alloc >= 0)
  expect_true(is.finite(res$kkt_boundary) && res$kkt_boundary >= 0)
  # the reported combined residual is the larger of the two
  expect_equal(res$kkt_residual, max(res$kkt_alloc, res$kkt_boundary))

  # per-start log is the right length and the best matches the objective
  expect_length(res$objective_by_start, 4)
  expect_length(res$converged_by_start, 4)
  expect_true(res$n_starts_failed >= 0L)
  ok <- res$objective_by_start[is.finite(res$objective_by_start)]
  expect_equal(res$objective, min(ok))
})

test_that("parametric mode runs and reports mode in .meta", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  dl  <- list(X1 = list(family = "lnorm", meanlog = 1.0, sdlog = 0.5),
              X2 = list(family = "lnorm", meanlog = 1.5, sdlog = 0.8))
  res <- stratify_multi(
    data = pop[, c("X1", "X2")], study_vars = pop[, c("Y1", "Y2")],
    H_vec = c(3, 3), n_total = 300, n_starts = 3, max_iter = 25,
    distributions = dl)
  expect_identical(res$.meta$mode, "parametric")
  expect_equal(sum(res$stratum_stats$W_h), 1, tolerance = 1e-8)
})

test_that("SLSQP method requires parametric mode and runs under it", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  expect_error(
    stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                   H_vec = c(3, 3), n_total = 300, method = "slsqp"),
    regexp = "parametric")

  dl  <- list(X1 = list(family = "lnorm", meanlog = 1.0, sdlog = 0.5),
              X2 = list(family = "lnorm", meanlog = 1.5, sdlog = 0.8))
  res <- stratify_multi(
    pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
    H_vec = c(3, 3), n_total = 300, n_starts = 2, max_iter = 15,
    method = "slsqp", distributions = dl)
  expect_s3_class(res, "multistrata")
  expect_equal(sum(res$allocation), 300)
})

test_that("analytic boundary gradient matches finite differences", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  dat <- as.matrix(pop[, c("X1", "X2")])
  sv  <- as.matrix(pop[, c("Y1", "Y2")])
  dl  <- list(X1 = list(family = "lnorm", meanlog = 1.0, sdlog = 0.5),
              X2 = list(family = "lnorm", meanlog = 1.5, sdlog = 0.8))
  H_vec <- c(3L, 3L)
  alpha <- apply(dat, 2, min);  beta <- apply(dat, 2, max)
  B     <- list(quantile(dat[, 1], c(1/3, 2/3), names = FALSE),
                quantile(dat[, 2], c(1/3, 2/3), names = FALSE))
  n_al  <- rep(300 / 9, 9)
  omega <- c(0.5, 0.5)

  for (j in 1:2) {
    ga <- multistrata:::.phi_grad_boundary(
      B[[j]], j, B, dat, sv, H_vec, alpha, beta, omega, n_al, dl, nrow(dat))
    # FD of the *smooth part*: parametric weights with sigma2 frozen at B
    sv0 <- multistrata:::.compute_stratum_means_vars(dat, sv, B, H_vec,
                                                        dl, alpha, beta)
    phi_W <- function(bj) {
      Bt <- B;  Bt[[j]] <- sort(bj)
      ps <- compute_stratum_stats_parametric(dat, NULL, Bt, H_vec,
                                             alpha, beta, dl)
      sum(omega * colSums(sv0$sigma2 * ps$W_h^2 *
                          (1 / n_al - 1 / pmax(ps$N_h, 1))))
    }
    gf <- ga * 0
    for (l in seq_along(B[[j]])) {
      e <- 1e-5 * (beta[j] - alpha[j])
      bp <- B[[j]]; bp[l] <- bp[l] + e
      bm <- B[[j]]; bm[l] <- bm[l] - e
      gf[l] <- (phi_W(bp) - phi_W(bm)) / (2 * e)
    }
    # N_h rounding in compute_stratum_stats_parametric adds small noise;
    # require close agreement relative to gradient scale
    expect_lt(max(abs(ga - gf)) / max(abs(gf), 1e-12), 0.05)
  }
})

test_that("CV-constrained joint design achieves attainable targets", {
  skip_if_not_installed("nloptr")
  pop  <- make_pop()
  base <- stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                         H_vec = c(3, 3), n_total = 400, n_starts = 2,
                         max_iter = 15)
  ss   <- base$stratum_stats
  s2   <- as.matrix(ss[, c("sigma2_k1", "sigma2_k2")])
  mu   <- as.matrix(ss[, c("mu_k1", "mu_k2")])
  cv0  <- design_cv(s2, mu, ss$W_h, base$allocation, ss$N_h)

  res <- stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                        H_vec = c(3, 3), n_total = 400, n_starts = 2,
                        max_iter = 15,
                        cv_targets = as.numeric(cv0) * 1.05)
  expect_false(is.null(res$cv_achieved))
  expect_true(all(res$cv_achieved <= as.numeric(cv0) * 1.05 + 5e-4))
})

test_that("mixed-radix stratum encoding matches expand.grid label order", {
  H_vec  <- c(3L, 4L)
  combos <- expand.grid(seq_len(3), seq_len(4))
  mult   <- cumprod(c(1L, H_vec[-length(H_vec)]))
  idx    <- 1L + as.integer(as.matrix(combos - 1L) %*% mult)
  expect_identical(idx, seq_len(12L))
})
