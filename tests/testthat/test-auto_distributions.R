# Tests for the auto-parametric mode (fit_distributions) --------------------

test_that("fit_distributions recovers the true families by AIC", {
  skip_on_cran()
  set.seed(101)
  X <- data.frame(
    A = rlnorm(3000, meanlog = 1.0, sdlog = 0.75),
    B = rgamma(3000, shape = 2, rate = 1.5),
    C = runif(3000, 2, 7)
  )
  dl <- fit_distributions(X, engine = "builtin")
  ch <- attr(dl, "chosen")
  expect_identical(unname(ch["A"]), "lnorm")
  expect_identical(unname(ch["B"]), "gamma")
  expect_identical(unname(ch["C"]), "unif")
  # returned specs are directly usable by the parametric machinery
  for (j in seq_along(dl)) {
    Fj <- multistrata:::.build_cdf(dl[[j]], min(X[[j]]), max(X[[j]]))
    expect_true(all(diff(Fj(seq(min(X[[j]]), max(X[[j]]), length.out = 50))) >= 0))
  }
  ft <- attr(dl, "fit_table")
  expect_true(all(c("variable", "family", "AIC", "engine") %in% names(ft)))
})

test_that("Pareto tail data selects the pareto family", {
  skip_on_cran()
  set.seed(102)
  x <- data.frame(P = 1.5 * (1 - runif(4000))^(-1 / 2.8))
  ch <- attr(fit_distributions(x, engine = "builtin"), "chosen")
  expect_identical(unname(ch["P"]), "pareto")
})

test_that("stratify_multi(distributions = 'auto') runs end-to-end", {
  skip_on_cran()
  skip_if_not_installed("nloptr")
  set.seed(103)
  pop <- sim_population(N = 1500, p = 2, q = 2, scenario = "skewed")
  res <- stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
                        H_vec = c(3, 3), n_total = 300,
                        n_starts = 2, max_iter = 15,
                        distributions = "auto")
  expect_identical(res$.meta$mode, "parametric")
  expect_false(is.null(res$.meta$auto_fit))
  expect_length(res$.meta$distributions, 2)
  expect_equal(sum(res$allocation), 300)
})

test_that("fitverse engine errors cleanly when FitVerse is absent", {
  skip_on_cran()
  if (requireNamespace("FitVerse", quietly = TRUE)) {
    dl <- fit_distributions(
      data.frame(A = rlnorm(2000, 1, 0.6)), engine = "fitverse")
    expect_true(attr(dl, "chosen")[["A"]] %in%
                  c("lnorm", "gamma", "weibull"))
  } else {
    expect_error(
      fit_distributions(data.frame(A = rlnorm(500, 1, 0.6)),
                        engine = "fitverse"),
      regexp = "FitVerse")
  }
})

# ---------------------------------------------------------------------------
# RNG hygiene around the large-column subsample.
#
# fit_distributions() seeds a fixed subsample when a column is longer than
# max_fit_n, so that AIC family selection does not wobble between runs. It
# must put the caller's random stream back afterwards. Getting this wrong is
# not cosmetic: stratify_multi(distributions = "auto") calls this function
# before drawing its multi-start boundaries, so a leaked seed silently
# overrides the user's own set.seed() for the entire optimisation.
# ---------------------------------------------------------------------------

seed_fixture <- function(N = 12000) {
  set.seed(777)
  sim_population(N = N, p = 2, q = 2, scenario = "correlated")
}

test_that("fit_distributions leaves the caller's RNG stream untouched", {
  skip_on_cran()
  pop <- seed_fixture()
  X <- as.matrix(pop[, c("X1", "X2")])
  # N = 12000 > the default max_fit_n = 8000, so the subsample branch fires.
  expect_gt(nrow(X), 8000L)

  set.seed(99); reference <- rnorm(3)
  set.seed(99); invisible(fit_distributions(X)); after <- rnorm(3)
  expect_equal(after, reference, tolerance = 0)

  # Also when .Random.seed did not exist beforehand, the function must not
  # leave one behind.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  invisible(fit_distributions(X))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the caller's seed governs an auto-mode design", {
  skip_on_cran()
  pop <- seed_fixture()
  X <- as.matrix(pop[, c("X1", "X2")])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  run <- function(s) {
    set.seed(s)
    stratify_multi(X, Y, H_vec = c(2, 2), n_total = 200, omega = c(0.5, 0.5),
                   n_starts = 3, max_iter = 25, distributions = "auto",
                   verbose = FALSE)$objective
  }

  # Different seeds must give different multi-start outcomes. Before the RNG
  # fix these were byte-identical, because the internal set.seed() had
  # replaced the user's stream before any boundary was drawn.
  expect_false(isTRUE(all.equal(run(1), run(2))))
  # ... and the same seed must still reproduce exactly.
  expect_equal(run(1), run(1), tolerance = 0)
})

test_that("auto-mode family selection is stable across seeds", {
  skip_on_cran()
  pop <- seed_fixture()
  X <- as.matrix(pop[, c("X1", "X2")])
  fit <- function(s) { set.seed(s); fit_distributions(X) }

  f1 <- fit(1); f2 <- fit(2)
  # The point of the fixed subsample seed: two analysts with the same frame
  # and different seeds must not end up asserting different marginal
  # families. Selection and the fitted parameters are both deterministic.
  expect_identical(attr(f1, "chosen"), attr(f2, "chosen"))
  expect_equal(f1$X1, f2$X1, tolerance = 0)
  expect_equal(f1$X2, f2$X2, tolerance = 0)
})
