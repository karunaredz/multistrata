make_sweep_pop <- function(N = 1200) {
  set.seed(4242)
  sim_population(N = N, p = 2, q = 2, scenario = "correlated")
}

test_that("omega_sweep returns one row per grid point with the documented columns", {
  pop <- make_sweep_pop()
  g <- c(0.25, 0.5, 0.75)
  sw <- omega_sweep(
    data       = pop[, c("X1", "X2")],
    study_vars = pop[, c("Y1", "Y2")],
    H_vec      = c(2, 2),
    n_total    = 150,
    grid       = g,
    n_starts   = 1,
    max_iter   = 20)

  expect_s3_class(sw, "data.frame")
  expect_identical(nrow(sw), length(g))
  expect_identical(names(sw), c("w", "phi", "cv1", "cv2", "gap"))
  expect_true(all(vapply(sw, function(x) all(is.finite(x)), logical(1))))
  expect_true(all(sw$phi > 0))
  expect_true(all(sw$cv1 > 0 & sw$cv2 > 0))
  expect_true(all(sw$gap >= 0))

  # Returned columns match omega_sweep_density()'s, so the two are
  # interchangeable downstream.
  expect_identical(attr(sw, "pair"), c(1L, 2L))
  expect_identical(attr(sw, "pair_names"), c("Y1", "Y2"))
})

test_that("omega_sweep renormalises weights so w tracks the requested share", {
  pop <- make_sweep_pop()
  g <- c(0.2, 0.8)
  sw <- omega_sweep(
    data       = pop[, c("X1", "X2")],
    study_vars = pop[, c("Y1", "Y2")],
    H_vec      = c(2, 2),
    n_total    = 150,
    omega      = c(0.5, 0.5),
    grid       = g,
    n_starts   = 1,
    max_iter   = 20)

  # With two study variables the pair carries all the weight, so the
  # renormalised weight is exactly the requested share.
  expect_equal(sw$w, g, tolerance = 1e-8)
})

test_that("omega_sweep trades precision between the two study variables", {
  set.seed(23)
  pop <- make_sweep_pop()
  sw <- omega_sweep(
    data       = pop[, c("X1", "X2")],
    study_vars = pop[, c("Y1", "Y2")],
    H_vec      = c(2, 2),
    n_total    = 200,
    grid       = c(0.1, 0.9),
    n_starts   = 2,
    max_iter   = 30)

  # Shifting priority onto Y1 must not make Y1's CV worse while also
  # improving Y2's: the sweep exists to expose a trade-off.
  expect_lte(sw$cv1[2], sw$cv1[1])
  expect_gte(sw$cv2[2], sw$cv2[1])
})

test_that("omega_sweep agrees with a hand-rolled stratify_multi sweep", {
  pop <- make_sweep_pop()
  X <- as.matrix(pop[, c("X1", "X2")])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  g <- c(0.3, 0.7)

  # The multi-start scheme draws random initial boundaries, so both the sweep
  # and the manual loop must start from the same RNG state for the comparison
  # to isolate omega_sweep()'s bookkeeping from the restart lottery.
  set.seed(11)
  sw <- omega_sweep(data = X, study_vars = Y, H_vec = c(2, 2), n_total = 150,
                    grid = g, n_starts = 1, max_iter = 20)

  set.seed(11)
  manual <- vapply(g, function(w) {
    r <- stratify_multi(data = X, study_vars = Y, H_vec = c(2, 2),
                        n_total = 150, omega = c(w, 1 - w),
                        n_starts = 1, max_iter = 20, verbose = FALSE)
    r$objective
  }, numeric(1))

  # omega_sweep() must be exactly a loop over stratify_multi(), not an
  # approximation of one.
  expect_equal(sw$phi, manual, tolerance = 1e-12)
})

test_that("omega_sweep accepts study-variable names and a fitted design", {
  pop <- make_sweep_pop()
  X <- as.matrix(pop[, c("X1", "X2")])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  fit <- stratify_multi(data = X, study_vars = Y, H_vec = c(2, 2),
                        n_total = 150, omega = c(0.5, 0.5),
                        n_starts = 1, max_iter = 20, verbose = FALSE)

  set.seed(7)
  by_name <- omega_sweep(data = X, study_vars = Y, n_total = 150,
                         pair = c("Y1", "Y2"), grid = 0.5,
                         result = fit, n_starts = 1, max_iter = 20)
  set.seed(7)
  by_index <- omega_sweep(data = X, study_vars = Y, n_total = 150,
                          pair = c(1L, 2L), grid = 0.5,
                          result = fit, n_starts = 1, max_iter = 20)

  # Naming the pair and indexing it must resolve to the same sweep.
  expect_equal(by_name$phi, by_index$phi, tolerance = 1e-12)
  expect_identical(attr(by_name, "pair_names"), c("Y1", "Y2"))

  # H_vec, omega and method are inherited from the fitted object.
  expect_identical(nrow(by_name), 1L)
})

test_that("omega_sweep validates its pair and grid arguments", {
  pop <- make_sweep_pop(400)
  X <- as.matrix(pop[, c("X1", "X2")])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  args <- list(data = X, study_vars = Y, H_vec = c(2, 2), n_total = 100,
               n_starts = 1, max_iter = 5)

  expect_error(do.call(omega_sweep, c(args, list(pair = c(1L, 1L)))),
               "two different study variables")
  expect_error(do.call(omega_sweep, c(args, list(pair = c(1L, 9L)))),
               "two different study variables")
  expect_error(do.call(omega_sweep, c(args, list(pair = c("Y1", "nope")))),
               "not found among the study variables")
  expect_error(do.call(omega_sweep, c(args, list(grid = c(0, 0.5)))),
               "strictly between 0 and 1")
  expect_error(do.call(omega_sweep, c(args, list(grid = 1.5))),
               "strictly between 0 and 1")
  expect_error(do.call(omega_sweep, c(args, list(omega = c(1, 1, 1)))),
               "must have length equal to")
  expect_error(omega_sweep(data = X, study_vars = Y, n_total = 100,
                           grid = 0.5, n_starts = 1),
               "'H_vec' must be supplied")
  expect_error(omega_sweep(data = X, study_vars = Y, H_vec = c(2, 2),
                           n_total = 100, grid = 0.5, result = "not a fit"),
               "must be a 'multistrata' object")
})

test_that("omega_sweep_density reports the swept pair the same way", {
  dists <- list(X1 = list(family = "exp", rate = 1),
                X2 = list(family = "lnorm", meanlog = 0, sdlog = 0.5))
  sw <- omega_sweep_density(
    distributions = dists, H_vec = c(2, 2), n_total = 200, N = 2000,
    grid = c(0.3, 0.7), n_starts = 1, max_iter = 20)

  expect_identical(names(sw), c("w", "phi", "cv1", "cv2", "gap"))
  expect_identical(attr(sw, "pair"), c(1L, 2L))
  expect_identical(attr(sw, "pair_names"), c("X1", "X2"))
})