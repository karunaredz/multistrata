# Tests for composite variance, CVs and the Cauchy-Schwarz lower bound ------

test_that("composite_variance matches a hand computation", {
  sigma2 <- matrix(c(4, 9, 1, 4), nrow = 2)
  W_h    <- c(0.4, 0.6)
  n_h    <- c(80, 120)
  N_h    <- c(400, 600)
  omega  <- c(0.5, 0.5)
  fpc    <- 1 / n_h - 1 / N_h
  manual <- 0.5 * sum(W_h^2 * sigma2[, 1] * fpc) +
            0.5 * sum(W_h^2 * sigma2[, 2] * fpc)
  expect_equal(composite_variance(sigma2, W_h, n_h, N_h, omega), manual)
})

test_that("lower bound is tight exactly at the (continuous) Neyman allocation", {
  set.seed(21)
  H <- 10
  sigma2 <- matrix(runif(H * 2, 1, 25), H, 2)
  W_h    <- as.numeric(prop.table(runif(H) + 0.2))
  N_h    <- rep(1e7, H)                        # FPC negligible, bounds off
  omega  <- c(0.3, 0.7)
  n_tot  <- 700

  A_h <- multistrata:::.compute_A_h(sigma2, W_h, omega)
  n_c <- n_tot * A_h / sum(A_h)                # continuous Neyman
  phi <- composite_variance(sigma2, W_h, n_c, N_h, omega)
  lb  <- variance_lower_bound(sigma2, W_h, N_h, omega, n_tot)
  expect_equal(phi, lb, tolerance = 1e-8)

  # any other allocation is strictly above the bound
  n_bad <- rep(n_tot / H, H)
  expect_gt(composite_variance(sigma2, W_h, n_bad, N_h, omega), lb)
})

test_that("design_cv matches a hand computation", {
  sigma2 <- matrix(c(4, 9, 1, 4), nrow = 2)
  mu     <- matrix(c(10, 20, 5, 15), nrow = 2)
  W_h    <- c(0.4, 0.6)
  n_h    <- c(80, 120)
  N_h    <- c(400, 600)
  fpc    <- 1 / n_h - 1 / N_h
  V1     <- sum(W_h^2 * sigma2[, 1] * fpc)
  mu1    <- sum(W_h * mu[, 1])
  cvs    <- design_cv(sigma2, mu, W_h, n_h, N_h)
  expect_equal(as.numeric(cvs[1]), sqrt(V1) / mu1)
})
