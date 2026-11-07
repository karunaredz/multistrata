# Tests for parametric stratum weights, CDF/PDF/quantile builders -----------

test_that("parametric weights sum to one and match empirical proportions", {
  set.seed(11)
  pop <- sim_population(N = 8000, p = 2, q = 2, scenario = "skewed")
  B   <- list(c(2.0, 4.5), c(2.5, 5.5))

  st <- compute_stratum_stats_parametric(
    data       = pop[, c("X1", "X2")],
    study_vars = pop[, c("Y1", "Y2")],
    B          = B,
    H_vec      = c(3L, 3L),
    alpha      = apply(pop[, c("X1", "X2")], 2, min),
    beta       = apply(pop[, c("X1", "X2")], 2, max),
    dist_list  = list(
      X1 = list(family = "lnorm", meanlog = 1.0, sdlog = 0.5),
      X2 = list(family = "lnorm", meanlog = 1.5, sdlog = 0.8)
    )
  )
  expect_equal(sum(st$W_h), 1, tolerance = 1e-10)
  # correctly specified family => analytic and empirical weights agree
  emp <- st$counts / sum(st$counts)
  expect_lt(max(abs(st$W_h - emp)), 0.03)
})

test_that("quantile builders invert the CDF builders for all 13 families", {
  fams <- list(
    list(family = "norm",    mean = 2, sd = 1.5),
    list(family = "lnorm",   meanlog = 0.5, sdlog = 0.8),
    list(family = "gamma",   shape = 2, rate = 1.5),
    list(family = "exp",     rate = 1.3),
    list(family = "weibull", shape = 1.7, scale = 2.2),
    list(family = "beta",    shape1 = 2, shape2 = 3),
    list(family = "unif",    min = 1, max = 9),
    list(family = "chisq",   df = 4),
    list(family = "t",       df = 8),
    list(family = "cauchy",  location = 0, scale = 1),
    list(family = "logis",   location = 1, scale = 0.7),
    list(family = "pareto",  shape = 2.5, scale = 1.2),
    list(family = "triang",  a = 0, b = 2, c = 2)
  )
  u <- seq(0.05, 0.95, by = 0.05)
  for (d in fams) {
    Fj <- multistrata:::.build_cdf(d, 0, 10)
    Qj <- multistrata:::.build_quantile(d, 0, 10)
    expect_equal(Fj(Qj(u)), u, tolerance = 1e-8,
                 info = paste("family:", d$family))
  }
})

test_that("PDF builders match finite differences of the CDF builders", {
  fams <- list(
    list(family = "lnorm",  meanlog = 0.5, sdlog = 0.8),
    list(family = "gamma",  shape = 2, rate = 1.5),
    list(family = "exp",    rate = 1.0),
    list(family = "pareto", shape = 2.5, scale = 1.2),
    list(family = "triang", a = 0, b = 2, c = 2),
    list(family = "beta",   shape1 = 2, shape2 = 3)
  )
  for (d in fams) {
    Fj <- multistrata:::.build_cdf(d, 0, 4)
    fj <- multistrata:::.build_pdf(d, 0, 4)
    xs <- seq(0.3, 3.5, by = 0.2)
    if (d$family == "pareto") xs <- xs[xs > 1.2]
    if (d$family == "triang") xs <- xs[xs < 2]
    fd <- (Fj(xs + 1e-6) - Fj(xs - 1e-6)) / 2e-6
    expect_equal(fj(xs), fd, tolerance = 1e-4,
                 info = paste("family:", d$family))
  }
})

test_that("right-triangular CDF matches the closed form of the paper", {
  # Paper eq. (27): decreasing right-triangular, density 2(b - x)/(b - a)^2 with
  # mode at a = 0, so F(x) = 1 - ((b - x)/(b - a))^2.  Mode at a means c = 0.
  # (c = 2 would be the mirror image, F(x) = (x/2)^2, and is covered by the
  #  builder round-trip tests above.)
  Fj <- multistrata:::.build_cdf(list(family = "triang", a = 0, b = 2, c = 0), 0, 2)
  xs <- seq(0, 2, by = 0.1)
  expect_equal(Fj(xs), 1 - ((2 - xs) / 2)^2, tolerance = 1e-12)
})
