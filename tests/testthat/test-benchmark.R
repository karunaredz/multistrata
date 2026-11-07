bench_fixture <- function(p = 2, H = c(2, 2), n_total = 200, N = 1500) {
  set.seed(808)
  pop <- sim_population(N = N, p = p, q = 2, scenario = "correlated")
  X <- as.matrix(pop[, paste0("X", seq_len(p)), drop = FALSE])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  res <- stratify_multi(X, Y, H_vec = H, n_total = n_total,
                        omega = c(0.5, 0.5), n_starts = 3, max_iter = 40,
                        verbose = FALSE)
  list(X = X, Y = Y, res = res, H = H, n_total = n_total)
}

test_that("benchmark_designs reports all four designs with consistent columns", {
  f <- bench_fixture()
  bm <- benchmark_designs(f$X, f$Y, H_vec = f$H, n_total = f$n_total,
                          omega = c(0.5, 0.5), result_joint = f$res)

  expect_s3_class(bm, "data.frame")
  expect_identical(names(bm), c("method", "composite_variance",
                                "relative_efficiency", "pct_gain_over_joint"))
  expect_identical(bm$method,
                   c("Joint (multistrata)", "Proportional allocation",
                     "Equal allocation", "Sequential (quantile + Neyman)"))
  expect_true(all(is.finite(bm$composite_variance)))
  expect_true(all(bm$composite_variance > 0))

  # relative_efficiency is phi_method / phi_joint, so the joint row is 1 and
  # pct_gain_over_joint is its percentage restatement.
  expect_equal(bm$relative_efficiency[1], 1, tolerance = 1e-8)
  expect_equal(bm$pct_gain_over_joint,
               round(100 * (bm$relative_efficiency - 1), 2),
               tolerance = 1e-8)
  expect_equal(bm$relative_efficiency,
               round(bm$composite_variance / bm$composite_variance[1], 4),
               tolerance = 1e-3)
})

test_that("the joint design is no worse than the benchmark alternatives", {
  f <- bench_fixture()
  bm <- benchmark_designs(f$X, f$Y, H_vec = f$H, n_total = f$n_total,
                          omega = c(0.5, 0.5), result_joint = f$res)

  # This is the package's operational claim, and Lemma 1 of the methodology:
  # solving boundaries and allocation jointly cannot be beaten by fixing
  # boundaries first (sequential) or by ignoring the variances
  # (proportional / equal allocation).
  phi_joint <- bm$composite_variance[1]
  expect_true(all(bm$composite_variance[-1] >= phi_joint - 1e-10))
  expect_true(all(bm$relative_efficiency[-1] >= 1 - 1e-6))
})

test_that("benchmark_designs runs the joint method itself when not supplied", {
  f <- bench_fixture()
  expect_message(
    bm <- benchmark_designs(f$X, f$Y, H_vec = f$H, n_total = f$n_total,
                            omega = c(0.5, 0.5), n_starts = 2),
    "Running joint stratify_multi")
  expect_identical(nrow(bm), 4L)
  expect_equal(bm$relative_efficiency[1], 1, tolerance = 1e-8)
})

test_that("relative_efficiency is the ratio of composite variances", {
  expect_equal(relative_efficiency(2, 4), 2)
  expect_equal(relative_efficiency(4, 4), 1)
  expect_equal(relative_efficiency(0.5, 0.25), 0.5)
  expect_error(relative_efficiency(0, 1), "must be positive")
  expect_error(relative_efficiency(-1, 1), "must be positive")
})

test_that("the optimal design respects the Cauchy-Schwarz lower bound", {
  f <- bench_fixture()
  ss <- f$res$stratum_stats
  q  <- 2L
  sigma2 <- as.matrix(ss[, paste0("sigma2_k", seq_len(q)), drop = FALSE])

  lb <- variance_lower_bound(sigma2, ss$W_h, ss$N_h, c(0.5, 0.5), f$n_total)
  expect_true(is.finite(lb))
  expect_gte(lb, 0)

  # The attained objective cannot fall below the analytic bound, and the
  # package's reported optimality gap must be that comparison.
  expect_gte(f$res$objective, lb - 1e-10)
  expect_equal(f$res$optimality_gap, (f$res$objective - lb) / lb,
               tolerance = 1e-6)
})