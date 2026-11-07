test_that("density benchmark compares analytic alternatives", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )
  r <- stratify_multi_density(
    d, H_vec = c(2, 2), n_total = 80, N = 500,
    ranges = list(c(0, 6), c(0, 2)),
    n_starts = 1, max_iter = 6
  )

  b <- benchmark_density_designs(result_joint = r)

  expect_s3_class(b, "data.frame")
  expect_equal(nrow(b), 4L)
  expect_true(all(c("method", "composite_variance", "relative_efficiency") %in% names(b)))
  expect_equal(b$relative_efficiency[1], 1)
  expect_true(all(is.finite(b$composite_variance)))
})

test_that("find_n_density locates a feasible sample size", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )
  r <- stratify_multi_density(
    d, H_vec = c(2, 2), n_total = 80, N = 500,
    ranges = list(c(0, 6), c(0, 2)),
    n_starts = 1, max_iter = 6
  )

  fn <- find_n_density(result = r, cv_targets = c(0.2, 0.2))

  expect_s3_class(fn, "findN")
  expect_true(fn$feasible)
  expect_gte(fn$n, prod(r$.meta$H_vec) * 2)
  expect_lte(fn$n, r$.meta$N)
  expect_true(all(fn$cv_achieved <= fn$cv_targets + 5e-4))
})

test_that("choose_H_density returns the same shape as choose_H", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )

  ch <- choose_H_density(
    d, n_total = 80, N = 500, H_each = 2:3,
    ranges = list(c(0, 6), c(0, 2)),
    n_starts = 1, max_iter = 5, plot = FALSE
  )

  expect_s3_class(ch, "chooseH")
  expect_equal(nrow(ch$table), 2L)
  expect_true(all(c("H_label", "H_total", "objective", "marginal_gain") %in% names(ch$table)))
  expect_length(ch$recommended, 2L)
})

test_that("omega_sweep_density returns objective and CV diagnostics", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )

  sw <- omega_sweep_density(
    d, H_vec = c(2, 2), n_total = 80, N = 500,
    ranges = list(c(0, 6), c(0, 2)),
    grid = c(0.25, 0.5, 0.75),
    n_starts = 1, max_iter = 5
  )

  expect_s3_class(sw, "data.frame")
  expect_equal(nrow(sw), 3L)
  expect_true(all(c("w", "phi", "cv1", "cv2", "gap") %in% names(sw)))
  expect_true(all(is.finite(sw$phi)))
  expect_equal(sw$w, c(0.25, 0.5, 0.75), tolerance = 1e-8)
})

test_that("density helpers handle unattainable CV targets gracefully", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )
  r <- stratify_multi_density(
    d, H_vec = c(2, 2), n_total = 60, N = 300,
    ranges = list(c(0, 6), c(0, 2)),
    n_starts = 1, max_iter = 5
  )

  fn <- NULL
  expect_warning(
    fn <- find_n_density(result = r, cv_targets = c(1e-6, 1e-6), n_range = c(8, 20)),
    "unattainable"
  )

  expect_s3_class(fn, "findN")
  expect_false(fn$feasible)
  expect_true(all(c("n", "feasible") %in% names(fn$trace)))
})

test_that("density analysis supports cost-constrained designs", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "triang", a = 0, b = 2, c = 0)
  )
  r <- stratify_multi_density(
    d, H_vec = c(2, 2), n_total = 80, N = 500,
    ranges = list(c(0, 6), c(0, 2)),
    constraint = "cost", cost_total = 90, cost_fixed = 10,
    cost_per_unit = 1, n_starts = 1, max_iter = 5
  )

  b <- benchmark_density_designs(
    result_joint = r, constraint = "cost", cost_total = 90,
    cost_fixed = 10, cost_per_unit = 1
  )

  expect_s3_class(b, "data.frame")
  expect_equal(nrow(b), 4L)
  expect_true(all(is.finite(b$relative_efficiency)))
  expect_lte(sum(r$allocation) + 10, 90)
})

test_that("omega_sweep_density validates pairs and supports three variables", {
  d <- list(
    X1 = list(family = "exp", rate = 1),
    X2 = list(family = "gamma", shape = 2, rate = 1),
    X3 = list(family = "triang", a = 0, b = 2, c = 0)
  )

  expect_error(
    omega_sweep_density(
      d, H_vec = c(2, 2, 2), n_total = 120, N = 600,
      ranges = list(c(0, 6), c(0, 8), c(0, 2)),
      pair = c(1, 1), n_starts = 1, max_iter = 4
    ),
    "two different"
  )

  sw <- omega_sweep_density(
    d, H_vec = c(2, 2, 2), n_total = 120, N = 600,
    ranges = list(c(0, 6), c(0, 8), c(0, 2)),
    omega = c(0.4, 0.35, 0.25), pair = c(1, 3),
    grid = c(0.2, 0.8), n_starts = 1, max_iter = 4
  )

  expect_s3_class(sw, "data.frame")
  expect_equal(nrow(sw), 2L)
  expect_true(all(is.finite(sw$phi)))
})
