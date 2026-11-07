# Tests for the generalised Neyman allocation (Theorem 1(b)/(c)) ------------

test_that("q = 1 reduces to classical Neyman n_h ~ W_h * sigma_h", {
  H      <- 5
  W_h    <- c(0.30, 0.25, 0.20, 0.15, 0.10)
  sigma  <- c(2, 5, 9, 14, 20)
  N_h    <- rep(10000, H)                      # bounds inactive
  n_tot  <- 1000

  alloc <- neyman_alloc(matrix(sigma^2, ncol = 1), W_h, omega = 1,
                        n_total = n_tot, N_h = N_h, n_min = 2)
  expected <- n_tot * (W_h * sigma) / sum(W_h * sigma)
  expect_equal(sum(alloc), n_tot)
  expect_true(max(abs(alloc - expected)) <= 1)  # integer rounding only
})

test_that("A_h is W_h * sqrt(sum_k omega_k sigma2_hk), not the old sweep", {
  H <- 9; q <- 2
  set.seed(1)
  sigma2 <- matrix(runif(H * q, 1, 25), H, q)
  W_h    <- as.numeric(prop.table(runif(H)))
  omega  <- c(0.5, 0.5)
  A_ref  <- W_h * sqrt(as.numeric(sigma2 %*% omega))
  expect_equal(as.numeric(multistrata:::.compute_A_h(sigma2, W_h, omega)),
               pmax(A_ref, 1e-12), tolerance = 1e-12)
})

test_that("stationarity: A_h^2 / n_h^2 constant on strata inside bounds", {
  set.seed(2)
  H <- 12
  sigma2 <- matrix(runif(H * 2, 1, 30), H, 2)
  W_h    <- as.numeric(prop.table(runif(H)))
  N_h    <- sample(200:900, H)
  alloc  <- neyman_alloc(sigma2, W_h, omega = c(0.4, 0.6),
                         n_total = 600, N_h = N_h, n_min = 2)
  A_h  <- multistrata:::.compute_A_h(sigma2, W_h, c(0.4, 0.6))
  # use the continuous solution for the stationarity check
  n_c  <- multistrata:::.active_set_n(A_h, 600, pmin(2, N_h), N_h)
  free <- n_c > 2 + 1e-9 & n_c < N_h - 1e-9
  lam  <- (A_h[free] / n_c[free])^2
  expect_lt((max(lam) - min(lam)) / max(lam), 1e-8)
  expect_equal(sum(alloc), 600)
})

test_that("cost mode: budget respected, n_h ~ A_h / sqrt(c_h), c0 honoured", {
  set.seed(3)
  H <- 8
  sigma2 <- matrix(runif(H * 2, 1, 16), H, 2)
  W_h    <- as.numeric(prop.table(runif(H)))
  N_h    <- rep(100000, H)                     # bounds inactive
  c_h    <- runif(H, 0.5, 4)
  omega  <- c(0.5, 0.5)

  alloc <- neyman_alloc(sigma2, W_h, omega, N_h = N_h, n_min = 0,
                        constraint = "cost", cost_total = 900,
                        cost_fixed = 100, cost_per_unit = c_h)
  expect_lte(sum(c_h * alloc), 800 + max(c_h))  # integer rounding slack

  A_h <- multistrata:::.compute_A_h(sigma2, W_h, omega)
  n_c <- multistrata:::.active_set_cost(A_h, 800, c_h,
                                           rep(0, H), N_h)
  # stationarity: A_h^2 / (c_h * n_h^2) constant
  lam <- A_h^2 / (c_h * n_c^2)
  expect_lt((max(lam) - min(lam)) / max(lam), 1e-8)
  expect_equal(sum(c_h * n_c), 800, tolerance = 1e-6)
})

test_that("cost-optimal beats the old sqrt(A/c) rule", {
  set.seed(4)
  H <- 12
  A_h <- runif(H, 0.1, 10)
  c_h <- runif(H, 0.5, 8)
  budget <- 800
  n_new  <- multistrata:::.active_set_cost(A_h, budget, c_h,
                                              rep(1e-9, H), rep(1e9, H))
  n_old  <- budget * sqrt(A_h / c_h) / sum(sqrt(A_h * c_h))
  expect_lt(sum(A_h^2 / n_new), sum(A_h^2 / n_old))
})

test_that("CV-constrained allocation meets attainable targets", {
  skip_if_not_installed("nloptr")
  set.seed(5)
  H <- 8; q <- 2
  sigma2 <- matrix(runif(H * q, 1, 20), H, q)
  W_h    <- as.numeric(prop.table(runif(H) + 0.3))
  mu     <- matrix(runif(H * q, 5, 20), H, q)
  N_h    <- sample(300:900, H)
  omega  <- c(0.5, 0.5)
  n_tot  <- 400

  base   <- neyman_alloc(sigma2, W_h, omega, n_total = n_tot,
                         N_h = N_h, n_min = 2)
  cv0    <- design_cv(sigma2, mu, W_h, base, N_h)
  gamma  <- as.numeric(cv0) * c(0.995, 1.20)   # mildly tighten k = 1

  alloc <- neyman_alloc(sigma2, W_h, omega, n_total = n_tot, N_h = N_h,
                        n_min = 2, cv_targets = gamma, mu = mu)
  cv_ach <- attr(alloc, "cv_achieved")
  expect_equal(sum(alloc), n_tot)
  expect_true(all(cv_ach <= gamma + 5e-4))     # integer-rounding slack
})

test_that("infeasible CV targets raise a warning", {
  skip_if_not_installed("nloptr")
  set.seed(6)
  H <- 6
  sigma2 <- matrix(runif(H * 2, 5, 20), H, 2)
  W_h    <- as.numeric(prop.table(runif(H) + 0.3))
  mu     <- matrix(runif(H * 2, 5, 20), H, 2)
  N_h    <- rep(400, H)
  expect_warning(
    neyman_alloc(sigma2, W_h, c(0.5, 0.5), n_total = 60, N_h = N_h,
                 n_min = 2, cv_targets = c(1e-5, 1e-5), mu = mu),
    regexp = "not attainable")
})

# ---------------------------------------------------------------------------
# CV-constrained allocation with empty strata.
#
# A trial boundary set during the search can leave a cell with no units. Such
# a stratum has N_h = 0 and W_h = 0 and contributes nothing to the variance,
# but writing the finite-population correction as 1/n - 1/N_h evaluates to Inf
# there, and W_h^2 * sigma2 * Inf is 0 * Inf = NaN. That NaN reached nloptr as
# its warm-start objective and aborted the run with "objective in x0 returns
# NA" -- so cv_targets failed more often as n_starts grew, and was unusable at
# the package default of 50.
# ---------------------------------------------------------------------------

test_that("the CV-constrained allocation tolerates empty strata", {
  # Stratum 2 is empty.
  sigma2 <- rbind(c(4, 9), c(1, 1), c(2, 3))
  mu     <- rbind(c(10, 20), c(5, 5), c(8, 8))
  W_h    <- c(0.6, 0.0, 0.4)
  N_h    <- c(600, 0, 400)

  alloc <- neyman_alloc(sigma2, W_h, omega = c(0.5, 0.5), n_total = 100,
                        N_h = N_h, n_min = 2,
                        cv_targets = c(0.20, 0.20), mu = mu)

  expect_true(all(is.finite(alloc)))
  expect_identical(sum(alloc), 100)
  # The empty stratum can take no sample.
  expect_identical(alloc[2], 0)
  expect_true(all(is.finite(attr(alloc, "cv_achieved"))))
})

test_that("cv_targets survives a multi-start run at the package default", {
  skip_on_cran()          # 50 starts is too slow for CRAN's time budget
  data(pacific_agri, package = "multistrata")
  X <- as.matrix(pacific_agri[, c("land_area", "labour_units")])
  Y <- as.matrix(pacific_agri[, c("crop_yield", "livestock_value",
                                  "expenditure")])
  set.seed(1)
  res <- stratify_multi(X, Y, H_vec = c(3, 3), n_total = 900,
                        omega = c(0.4, 0.4, 0.2),
                        cv_targets = c(0.05, 0.05, 0.10))

  expect_s3_class(res, "multistrata")
  expect_identical(sum(res$allocation), 900)
  expect_true(all(res$cv_achieved <= c(0.05, 0.05, 0.10) + 5e-4))
})

test_that("cv_targets runs across many independent starts", {
  # The original failure was start-dependent, so cover a spread of seeds
  # rather than a single one.
  data(pacific_agri, package = "multistrata")
  X <- as.matrix(pacific_agri[, c("land_area", "labour_units")])
  Y <- as.matrix(pacific_agri[, c("crop_yield", "livestock_value",
                                  "expenditure")])
  for (s in 1:4) {
    set.seed(s)
    expect_no_error(
      stratify_multi(X, Y, H_vec = c(3, 3), n_total = 900,
                     omega = c(0.4, 0.4, 0.2), n_starts = 4, max_iter = 20,
                     cv_targets = c(0.10, 0.10, 0.20)))
  }
})
