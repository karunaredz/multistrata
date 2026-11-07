# Tests for the operational helpers added in 1.1-0 --------------------------
#   find_n(), choose_H(), assign_strata(), as_svydesign()
# and the new sim_population() scenarios (weibull_pareto, copula).

make_pop <- function(N = 1500) {
  set.seed(31)
  sim_population(N = N, p = 2, q = 2, scenario = "skewed")
}

fit_small <- function(pop) {
  stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                 H_vec = c(3, 3), n_total = 300, n_starts = 3, max_iter = 20)
}

test_that("new sim_population scenarios produce valid populations", {
  for (sc in c("weibull_pareto", "copula")) {
    pop <- sim_population(N = 800, p = 2, q = 2, scenario = sc, seed = 7)
    expect_equal(nrow(pop), 800)
    expect_equal(ncol(pop), 4)
    expect_true(all(is.finite(as.matrix(pop))))
    expect_true(all(pop$X1 > 0) && all(pop$X2 > 0))   # positive supports
  }
  # copula margins should be genuinely correlated
  pc <- sim_population(N = 3000, p = 2, q = 2, scenario = "copula",
                       rho = 0.6, seed = 11)
  expect_gt(cor(pc$X1, pc$X2), 0.3)
})

test_that("assign_strata scores a frame consistently", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  res <- fit_small(pop)

  f <- assign_strata(res, pop[, c("X1", "X2")])
  expect_s3_class(f, "factor")
  expect_length(f, nrow(pop))
  expect_identical(levels(f), res$stratum_stats$stratum)

  # integer form must agree with the internal assignment used to fit
  i <- assign_strata(res, pop[, c("X1", "X2")], labels = FALSE)
  internal <- multistrata:::.assign_strata(as.matrix(pop[, c("X1", "X2")]),
                                              res$boundaries, res$.meta$H_vec)
  expect_identical(as.integer(i), internal)

  # column-name matching: reordered columns still score identically
  rev_df <- pop[, c("X2", "X1")]
  i2 <- assign_strata(res, rev_df, labels = FALSE)
  expect_identical(as.integer(i2), internal)

  # wrong number of unnamed columns errors
  unnamed <- unname(as.matrix(pop[, c("X1", "X2", "Y1")]))
  expect_error(assign_strata(res, unnamed),
               regexp = "columns")
})

test_that("as_svydesign builds a coherent frame", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  res <- fit_small(pop)

  s    <- assign_strata(res, pop[, c("X1", "X2")], labels = FALSE)
  take <- unlist(lapply(seq_along(res$allocation), function(h) {
    idx <- which(as.integer(s) == h)
    if (!length(idx)) return(integer(0))
    sample(idx, min(length(idx), res$allocation[h]))
  }))
  frame <- as_svydesign(res, pop[take, ], frame_only = TRUE)

  expect_true(all(c(".stratum", ".N_h", ".n_h", ".weight") %in% names(frame)))
  expect_equal(nrow(frame), length(take))
  # design weights are positive and finite
  expect_true(all(is.finite(frame$.weight) & frame$.weight > 0))
  # weight = N_h / n_h by construction
  expect_equal(frame$.weight, frame$.N_h / frame$.n_h, tolerance = 1e-9)

  if (requireNamespace("survey", quietly = TRUE)) {
    des <- as_svydesign(res, pop[take, ])
    expect_s3_class(des, "survey.design")
  }
})

test_that("find_n returns the smallest feasible n and is monotone", {
  skip_if_not_installed("nloptr")
  pop <- make_pop(2000)
  res <- fit_small(pop)

  # A loose target should be attainable at a modest n
  fn <- find_n(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
               cv_targets = c(0.08, 0.08), result = res)
  expect_s3_class(fn, "findN")
  expect_true(fn$feasible)
  expect_true(is.finite(fn$n) && fn$n >= prod(res$.meta$H_vec) * 2)
  expect_true(all(fn$cv_achieved <= c(0.08, 0.08) + 1e-6))

  # n - 1 should be infeasible (minimality), when above the lower bracket
  if (fn$n > prod(res$.meta$H_vec) * 2) {
    sv  <- multistrata:::.compute_stratum_means_vars(
      as.matrix(pop[, c("X1", "X2")]), as.matrix(pop[, c("Y1", "Y2")]),
      res$boundaries, res$.meta$H_vec, NULL,
      res$.meta$alpha, res$.meta$beta)
    al  <- suppressWarnings(neyman_alloc(
      sv$sigma2, sv$W_h, res$.meta$omega, n_total = fn$n - 1L, N_h = sv$N_h,
      n_min = 2L, cv_targets = c(0.08, 0.08), mu = sv$mu))
    expect_false(isTRUE(attr(al, "cv_feasible")) &&
                 all(attr(al, "cv_achieved") <= c(0.08, 0.08) + 1e-6))
  }
})

test_that("find_n flags unattainable targets", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  res <- fit_small(pop)
  expect_warning(
    fn <- find_n(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                 cv_targets = c(1e-5, 1e-5), result = res,
                 n_range = c(50, 400)),
    regexp = "unattainable")
  expect_false(fn$feasible)
  expect_true(is.na(fn$n))
})

test_that("choose_H traces Phi* and recommends a grid", {
  skip_if_not_installed("nloptr")
  pop <- make_pop()
  ch  <- choose_H(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                  n_total = 400, H_each = 2:4, n_starts = 3, max_iter = 20,
                  plot = FALSE)
  expect_s3_class(ch, "chooseH")
  expect_equal(nrow(ch$table), 3)
  # Phi* is non-increasing as the grid refines (more strata never hurt)
  expect_true(all(diff(ch$table$objective) <= 1e-8))
  expect_length(ch$recommended, 2)
  expect_true(ch$recommended_row %in% seq_len(3))
})

test_that("choose_H preserves cost constraints for every candidate", {
  skip_if_not_installed("nloptr")
  pop <- make_pop(600)
  ch <- choose_H(
    pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
    n_total = 120, H_each = 2:3, n_starts = 1, max_iter = 5,
    constraint = "cost", cost_total = 130, cost_fixed = 10,
    cost_per_unit = 1, plot = FALSE)

  expect_true(all(vapply(ch$results,
    function(x) identical(x$.meta$constraint, "cost"), logical(1))))
  expect_true(all(vapply(ch$results,
    function(x) sum(x$allocation) <= 120, logical(1))))
})
