# Tests for the data-free density mode: stratify_multi_density() ------------

test_that("quadrature stratum moments match closed forms (uniform)", {
  # Uniform(0,1) split at 0.5: each half has mean 0.25/0.75 and var 1/48.
  m <- multistrata:::.density_one_margin(
    list(family = "unif", min = 0, max = 1),
    cuts = 0.5, range = c(0, 1), Hj = 2L)
  expect_equal(m$p,    c(0.5, 0.5),               tolerance = 1e-6)
  expect_equal(m$mean, c(0.25, 0.75),             tolerance = 1e-6)
  expect_equal(m$var,  c(1/48, 1/48),             tolerance = 1e-6)
})

test_that("quadrature matches the exponential closed form", {
  # Exp(1) on [0,10], interval [0, 0.3917]:
  #   p    = 1 - e^-b
  #   mean = ( 1 - (b+1)e^-b ) / p
  b  <- 0.3917
  m  <- multistrata:::.density_one_margin(
    list(family = "exp", rate = 1), cuts = b, range = c(0, 10), Hj = 2L)
  p1    <- 1 - exp(-b)
  mean1 <- (1 - (b + 1) * exp(-b)) / p1
  # p is renormalised by the total mass on [0,10] (~1 - e^-10)
  expect_equal(m$mean[1], mean1, tolerance = 1e-5)
  expect_equal(m$p[1] / (m$p[1] + m$p[2]),
               p1 / (1 - exp(-10)), tolerance = 1e-5)
})

test_that("assembled weights factorise and sum to one", {
  marg <- list(
    list(p = c(0.5, 0.3, 0.2), mean = c(1, 2, 3), var = c(0.1, 0.2, 0.3)),
    list(p = c(0.6, 0.4),      mean = c(5, 6),     var = c(0.4, 0.5)))
  a <- multistrata:::.assemble_from_margins(marg, H_vec = c(3L, 2L), N = 1000)
  expect_equal(sum(a$W_h), 1, tolerance = 1e-12)
  expect_equal(length(a$W_h), 6L)
  # W for (i=1,j=1) = 0.5 * 0.6 = 0.30
  expect_equal(a$W_h[1], 0.30, tolerance = 1e-12)
  expect_true(all(a$N_h >= 1))
})

test_that("stratify_multi_density returns a coherent data-free solution", {
  skip_if_not_installed("nloptr")
  res <- stratify_multi_density(
    distributions = list(x = list(family = "exp",    rate = 1),
                         y = list(family = "triang", a = 0, b = 2, c = 0)),
    H_vec   = c(3, 2), n_total = 300, N = 2000,
    ranges  = list(c(0, 10), c(0, 2)),
    n_starts = 3, max_iter = 15)

  expect_s3_class(res, "multistrata")
  expect_identical(res$.meta$mode, "density")
  expect_equal(sum(res$stratum_stats$W_h), 1, tolerance = 1e-8)
  expect_equal(sum(res$allocation), 300)
  expect_length(res$boundaries, 2)
  # boundaries strictly increasing and inside their ranges
  b1 <- res$boundaries[[1]]; b2 <- res$boundaries[[2]]
  expect_true(all(diff(b1) > 0) && all(b1 > 0 & b1 < 10))
  expect_true(all(b2 > 0 & b2 < 2))
  expect_true(is.finite(res$objective) && res$objective > 0)
  expect_true(is.finite(res$optimality_gap) && res$optimality_gap >= -1e-8)
})

test_that("density weights agree with the parametric data-mode weights", {
  # For the same boundaries, analytic W_h from .assemble_from_margins must equal
  # compute_stratum_stats_parametric()'s W_h (both are prod_j Delta F_j).
  dl <- list(x = list(family = "exp",    rate = 1),
             y = list(family = "triang", a = 0, b = 2, c = 0))
  B  <- list(c(0.5, 1.2), c(0.8))
  H  <- c(3L, 2L)
  marg <- multistrata:::.density_margins(dl, B, H, list(c(0, 10), c(0, 2)), 200L)
  a    <- multistrata:::.assemble_from_margins(marg, H, N = 5000)

  # emulate a tiny frame only to call the parametric weight routine on the SAME
  # boundaries; weights are analytic so they do not depend on the sample.
  set.seed(1)
  x <- pmin(rexp(5000, 1), 10); y <- 2 - 2 * sqrt(1 - runif(5000))
  ps <- compute_stratum_stats_parametric(
    data = cbind(x, y), study_vars = NULL, B = B, H_vec = H,
    alpha = c(0, 0), beta = c(10, 2), dist_list = dl)
  expect_equal(sort(a$W_h), sort(ps$W_h), tolerance = 1e-6)
})
