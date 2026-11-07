# =============================================================================
# sim_population.R, unchanged logic, unchanged public API
# Synthetic population generator for testing and benchmarking multistrata.
# =============================================================================


#' Simulate a Synthetic Multivariate Survey Population
#'
#' Generates a population data frame with \eqn{p} stratification variables
#' \eqn{X_1, \ldots, X_p} and \eqn{q} study variables \eqn{Y_1, \ldots, Y_q}.
#'
#' @param N        Population size.  Default \code{10000}.
#' @param p        Number of stratification variables.  Default \code{2}.
#' @param q        Number of study variables.  Default \code{2}.
#' @param scenario One of \code{"independent"}, \code{"correlated"},
#'   \code{"skewed"}, \code{"heavy_tail"}, \code{"mixed"},
#'   \code{"weibull_pareto"} (Scenario C: heavy-tailed establishment-survey
#'   archetype, Weibull \eqn{\times} Pareto margins), or \code{"copula"}
#'   (Scenario E: dependent stratification variables via a Gaussian copula
#'   with log-normal margins).
#' @param rho      Correlation used in the \code{"correlated"} and
#'   \code{"copula"} scenarios (for \code{"copula"} it is the Gaussian-copula
#'   correlation between the latent normals).  Default \code{0.7}.
#' @param seed     Integer seed.  Default \code{42}.
#'
#' @return A \code{data.frame} with \eqn{N} rows and \eqn{p + q} columns
#'   named \code{X1, ..., Xp, Y1, ..., Yq}.
#'
#' @examples
#' pop <- sim_population(N = 5000, p = 2, q = 2, scenario = "independent")
#' head(pop)
#'
#' @export
sim_population <- function(N        = 10000L,
                            p        = 2L,
                            q        = 2L,
                            scenario = c("independent", "correlated",
                                         "skewed", "heavy_tail", "mixed",
                                         "weibull_pareto", "copula"),
                            rho      = 0.7,
                            seed     = 42L) {

  scenario <- match.arg(scenario)
  set.seed(seed)

  X <- switch(scenario,

    "independent" = {
      mat <- matrix(rnorm(N * p, mean = 5, sd = 2), nrow = N, ncol = p)
      abs(mat) + 0.1
    },

    "correlated" = {
      Sigma_X      <- matrix(rho, nrow = p, ncol = p)
      diag(Sigma_X) <- 1
      mat <- MASS::mvrnorm(N, mu = rep(5, p), Sigma = 4 * Sigma_X)
      abs(mat) + 0.1
    },

    "skewed" = {
      mat         <- matrix(0, nrow = N, ncol = p)
      meanlog_vec <- seq(1.0, 1.5, length.out = p)
      sdlog_vec   <- seq(0.5, 0.8, length.out = p)
      for (j in seq_len(p))
        mat[, j] <- rlnorm(N, meanlog = meanlog_vec[j], sdlog = sdlog_vec[j])
      mat
    },

    "heavy_tail" = {
      alpha_vec <- seq(2.5, 3.5, length.out = p)
      mat       <- matrix(0, nrow = N, ncol = p)
      for (j in seq_len(p)) {
        u        <- runif(N)
        mat[, j] <- u^(-1 / alpha_vec[j])
      }
      for (j in seq_len(p)) {
        cap           <- quantile(mat[, j], 0.995)
        mat[mat[, j] > cap, j] <- cap
      }
      mat
    },

    "mixed" = {
      mat      <- matrix(0, nrow = N, ncol = p)
      mat[, 1] <- rlnorm(N, meanlog = 1.2, sdlog = 0.6)
      if (p >= 2) mat[, 2] <- abs(rnorm(N, mean = 5, sd = 2)) + 0.1
      if (p >= 3) for (j in 3:p) mat[, j] <- runif(N, min = 1, max = 10)
      mat
    },

    # Scenario C: heavy tails.  Column 1 Weibull(shape 1.2, scale 2),
    # column 2 Pareto(shape 2.5, scale 1); further columns alternate.
    # Sparse tail cells are exactly where empirical stratum counts are
    # noisiest and the density-based (parametric) mode gains the most.
    "weibull_pareto" = {
      mat <- matrix(0, nrow = N, ncol = p)
      for (j in seq_len(p)) {
        if (j %% 2L == 1L) {
          mat[, j] <- rweibull(N, shape = 1.2, scale = 2)
        } else {
          u        <- runif(N)
          mat[, j] <- 1 * (1 - u)^(-1 / 2.5)          # Pareto(shape 2.5, scale 1)
        }
      }
      for (j in seq_len(p)) {                          # trim extreme top 0.5%
        cap <- quantile(mat[, j], 0.995)
        mat[mat[, j] > cap, j] <- cap
      }
      mat
    },

    # Scenario E: dependent stratification variables.  Gaussian copula
    # (equicorrelation rho on the latent normals) with log-normal margins.
    # Tests robustness of the independence factorisation W_h = prod_j DeltaF_j.
    "copula" = {
      Sigma_Z       <- matrix(rho, nrow = p, ncol = p)
      diag(Sigma_Z) <- 1
      Z             <- MASS::mvrnorm(N, mu = rep(0, p), Sigma = Sigma_Z)
      U             <- pnorm(Z)                          # uniform margins
      mat           <- matrix(0, nrow = N, ncol = p)
      meanlog_vec   <- seq(0.8, 1.2, length.out = p)
      sdlog_vec     <- seq(0.6, 0.8, length.out = p)
      for (j in seq_len(p))
        mat[, j] <- qlnorm(U[, j], meanlog = meanlog_vec[j], sdlog = sdlog_vec[j])
      mat
    }
  )

  coef_mat <- matrix(runif(p * q, min = 0.3, max = 1.5), nrow = p, ncol = q)
  Y_mean   <- X %*% coef_mat
  noise_sd <- sqrt(0.5 + 0.1 * abs(Y_mean))
  Y        <- Y_mean + matrix(rnorm(N * q), nrow = N, ncol = q) * noise_sd
  Y        <- abs(Y) + 0.01

  colnames(X) <- paste0("X", seq_len(p))
  colnames(Y) <- paste0("Y", seq_len(q))
  as.data.frame(cbind(X, Y))
}
