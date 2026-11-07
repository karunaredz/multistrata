# =============================================================================
# stratum_integrals.R
#
# Parametric density-based computation of stratum quantities:
#   W_h, stratum weight (mass under marginal density)
#   mu_hk, conditional mean of study variable Y_k in stratum h
#   sigma2_hk, conditional variance of Y_k in stratum h
#
# When the marginal distribution of each stratification variable X_j is
# specified parametrically, stratum weights factorise as
#
#   W_h = prod_j [ F_j(b_hj) - F_j(b_{h-1,j}) ]
#
# (Equation 2 in Reddy & Khan, forthcoming) which requires only marginal CDFs.
# This delivers a smoother objective landscape, removes dependence on finite
# population noise in W_h, and allows stratification without the full frame,
# only distributional parameters are needed.
#
# Study-variable moments (mu_hk, sigma2_hk) are still estimated empirically
# from the data, because a parametric model for Y | X is rarely available.
# The user may supply study_vars = NULL to work in a boundaries-only mode
# (no Y moments; only stratum weights and X-side statistics are returned).
#
# Supported distribution families (matched to base R d/p/q functions):
#   "norm"     Normal(mean, sd)
#   "lnorm"    Log-normal(meanlog, sdlog)
#   "gamma"    Gamma(shape, rate)  [or scale = 1/rate]
#   "exp"      Exponential(rate)
#   "weibull"  Weibull(shape, scale)
#   "beta"     Beta(shape1, shape2)  [support rescaled to (alpha_j, beta_j)]
#   "unif"     Uniform(min, max)
#   "chisq"    Chi-squared(df)
#   "t"        Student-t(df, ncp)
#   "cauchy"   Cauchy(location, scale)
#   "logis"    Logistic(location, scale)
#   "pareto"   Pareto(shape, scale), custom, not in base R stats
#   "triang"   Triangular(a, b, c), custom (right-triangular: c=b)
#
# For "beta" the support is linearly rescaled from (0,1) to (alpha_j, beta_j).
# For "pareto" and "triang" lightweight CDF/PDF helpers are provided below.
# =============================================================================


# -----------------------------------------------------------------------------
# Public: compute_stratum_stats_parametric
# -----------------------------------------------------------------------------

#' Parametric Stratum Statistics for Multivariate Stratification
#'
#' Computes stratum weights \eqn{W_h}, stratum population sizes \eqn{N_h},
#' stratum X-variable means and variances, and (empirically) study-variable
#' means and variances, using parametric marginal distributions for the
#' stratification variables.
#'
#' Stratum weights are computed analytically via the marginal CDFs (they do
#' not depend on the observed sample, only on the boundary configuration and
#' the specified distributions).  Study-variable statistics are estimated
#' from the data because a parametric model for Y|X is not assumed.
#'
#' @param data        Matrix/data.frame of stratification variables (N x p).
#' @param study_vars  Matrix/data.frame of study variables (N x q), or
#'   \code{NULL} to skip study-variable statistics.
#' @param B           List of length p; each element is a numeric vector of
#'   interior boundary points for the corresponding variable.
#' @param H_vec       Integer vector of length p (number of strata per variable).
#' @param alpha       Lower bound vector (length p); typically
#'   \code{apply(data, 2, min)}.
#' @param beta        Upper bound vector (length p); typically
#'   \code{apply(data, 2, max)}.
#' @param dist_list   Named list of length p.  Each element must be a list
#'   with a \code{family} character field and additional distribution
#'   parameters.  See Details.
#' @param n_alloc     Optional allocation vector (length H_total) used only to
#'   attach to the returned data.frame.  Not used in computations.
#'
#' @details
#' Each element of \code{dist_list} must contain at least \code{family}.
#' Additional fields match the base-R parameterisation:
#' \describe{
#'   \item{\code{"norm"}}{  \code{mean}, \code{sd}}
#'   \item{\code{"lnorm"}}{\code{meanlog}, \code{sdlog}}
#'   \item{\code{"gamma"}}{\code{shape}, \code{rate} (or \code{scale})}
#'   \item{\code{"exp"}}{  \code{rate}}
#'   \item{\code{"weibull"}}{\code{shape}, \code{scale}}
#'   \item{\code{"beta"}}{  \code{shape1}, \code{shape2}, rescaled to
#'     (\code{alpha_j}, \code{beta_j})}
#'   \item{\code{"unif"}}{ \code{min}, \code{max} (defaults to
#'     (\code{alpha_j}, \code{beta_j}) if omitted)}
#'   \item{\code{"chisq"}}{\code{df}}
#'   \item{\code{"t"}}{    \code{df}, optional \code{ncp}}
#'   \item{\code{"cauchy"}}{\code{location}, \code{scale}}
#'   \item{\code{"logis"}}{\code{location}, \code{scale}}
#'   \item{\code{"pareto"}}{\code{shape}, \code{scale}}
#'   \item{\code{"triang"}}{\code{a}, \code{b}, \code{c} (mode); set
#'     \code{c = b} for a right-triangular distribution)}
#' }
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{W_h}}{Numeric vector of stratum weights (analytically derived).}
#'   \item{\code{N_h}}{Effective population sizes \eqn{N_h = W_h \times N}.}
#'   \item{\code{mu}}{Matrix (H x q) of empirical study-variable stratum means.}
#'   \item{\code{sigma2}}{Matrix (H x q) of empirical study-variable variances.}
#'   \item{\code{strat_id}}{Integer vector of stratum assignments (length N).}
#'   \item{\code{counts}}{Integer vector of sample counts per stratum.}
#'   \item{\code{mode}}{Character "parametric".}
#' }
#'
#' @examples
#' set.seed(1)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "independent")
#' B   <- list(c(1.5, 3.0), c(0.5, 1.0))
#' stats <- compute_stratum_stats_parametric(
#'   data       = pop[, c("X1","X2")],
#'   study_vars = pop[, c("Y1","Y2")],
#'   B          = B,
#'   H_vec      = c(3L, 3L),
#'   alpha      = apply(pop[,c("X1","X2")], 2, min),
#'   beta       = apply(pop[,c("X1","X2")], 2, max),
#'   dist_list  = list(
#'     X1 = list(family = "lnorm", meanlog = 0.5, sdlog = 0.9),
#'     X2 = list(family = "exp",   rate    = 1.0)
#'   )
#' )
#' round(stats$W_h, 4)
#'
#' @export
compute_stratum_stats_parametric <- function(data,
                                              study_vars,
                                              B,
                                              H_vec,
                                              alpha,
                                              beta,
                                              dist_list,
                                              n_alloc = NULL) {

  data  <- as.matrix(data)
  N     <- nrow(data)
  p     <- ncol(data)
  q     <- if (!is.null(study_vars)) ncol(as.matrix(study_vars)) else 0L

  if (length(dist_list) != p)
    stop("'dist_list' must have one entry per column of 'data' (length p = ", p, ").")
  if (length(H_vec) != p)
    stop("'H_vec' must have length p = ", p, ".")

  # ------------------------------------------------------------------
  # 1.  Validate and pre-build CDF callables for each X_j
  # ------------------------------------------------------------------
  cdf_fns <- vector("list", p)
  for (j in seq_len(p)) {
    d <- dist_list[[j]]
    if (is.null(d$family))
      stop("dist_list[[", j, "]] must contain a 'family' field.")
    cdf_fns[[j]] <- .build_cdf(d, alpha[j], beta[j])
  }

  # ------------------------------------------------------------------
  # 2.  Enumerate all H_total strata and compute W_h analytically
  # ------------------------------------------------------------------
  H_total <- prod(H_vec)
  labels  <- .stratum_index_labels(H_vec)

  # Full boundary sequence per variable: (alpha_j, b1, ..., b_{Hj-1}, beta_j)
  full_B <- vector("list", p)
  for (j in seq_len(p))
    full_B[[j]] <- c(alpha[j], B[[j]], beta[j])

  # Enumerate multi-index h = (h1,...,hp) for all strata
  grids  <- lapply(seq_len(p), function(j) seq_len(H_vec[j]))
  combos <- expand.grid(grids, KEEP.OUT.ATTRS = FALSE)  # H_total x p

  W_h_param <- numeric(H_total)
  for (h in seq_len(H_total)) {
    wt <- 1.0
    for (j in seq_len(p)) {
      hj  <- combos[h, j]
      lo  <- full_B[[j]][hj]
      hi  <- full_B[[j]][hj + 1L]
      wt  <- wt * (cdf_fns[[j]](hi) - cdf_fns[[j]](lo))
    }
    W_h_param[h] <- max(wt, 0)
  }

  # Normalise to sum to 1 (numerical safety)
  s <- sum(W_h_param)
  if (s < 1e-12)
    stop("All parametric stratum weights are essentially zero. ",
         "Check that dist_list parameters are consistent with the data range.")
  W_h_param <- W_h_param / s

  N_h_param <- pmax(round(W_h_param * N), 1L)

  # ------------------------------------------------------------------
  # 3.  Empirical stratum assignments (for study-variable moments)
  # ------------------------------------------------------------------
  strat_id <- .assign_strata(data, B, H_vec)
  counts   <- tabulate(strat_id, nbins = H_total)

  # ------------------------------------------------------------------
  # 4.  Empirical study-variable moments (Y within strata)
  # ------------------------------------------------------------------
  if (q > 0L) {
    sv_mat <- as.matrix(study_vars)
    mu     <- matrix(0, nrow = H_total, ncol = q)
    sigma2 <- matrix(0, nrow = H_total, ncol = q)
    global_mu     <- colMeans(sv_mat)
    global_sigma2 <- apply(sv_mat, 2L, var)

    for (h in seq_len(H_total)) {
      idx <- which(strat_id == h)
      if (length(idx) < 2L) {
        mu[h, ]     <- global_mu
        sigma2[h, ] <- global_sigma2
      } else {
        ys          <- sv_mat[idx, , drop = FALSE]
        mu[h, ]     <- colMeans(ys)
        sigma2[h, ] <- apply(ys, 2L, var)
      }
    }
  } else {
    mu     <- matrix(NA_real_, nrow = H_total, ncol = 0L)
    sigma2 <- matrix(NA_real_, nrow = H_total, ncol = 0L)
  }

  list(
    W_h      = W_h_param,
    N_h      = N_h_param,
    mu       = mu,
    sigma2   = sigma2,
    strat_id = strat_id,
    counts   = counts,
    mode     = "parametric"
  )
}


# =============================================================================
# Internal: distribution helpers
# =============================================================================

#' Build a CDF function for a given distribution specification
#' @keywords internal
.build_cdf <- function(d, alpha_j, beta_j) {
  fam <- tolower(trimws(d$family))

  switch(fam,

    "norm" = {
      mn <- d$mean %||% 0;  sd <- d$sd %||% 1
      function(x) pnorm(x, mean = mn, sd = sd)
    },

    "lnorm" = {
      ml <- d$meanlog %||% 0;  sl <- d$sdlog %||% 1
      function(x) plnorm(x, meanlog = ml, sdlog = sl)
    },

    "gamma" = {
      sh <- d$shape
      if (is.null(sh)) stop("gamma distribution requires 'shape'.")
      rt <- if (!is.null(d$rate)) d$rate else if (!is.null(d$scale)) 1 / d$scale else 1
      function(x) pgamma(x, shape = sh, rate = rt)
    },

    "exp" = {
      rt <- d$rate %||% 1
      function(x) pexp(x, rate = rt)
    },

    "weibull" = {
      sh <- d$shape %||% 1;  sc <- d$scale %||% 1
      function(x) pweibull(x, shape = sh, scale = sc)
    },

    "beta" = {
      # Beta is on (0,1); linearly rescale to (alpha_j, beta_j)
      s1 <- d$shape1 %||% 1;  s2 <- d$shape2 %||% 1
      rng <- beta_j - alpha_j
      function(x) pbeta((x - alpha_j) / rng, shape1 = s1, shape2 = s2)
    },

    "unif" = {
      lo <- if (!is.null(d$min)) d$min else alpha_j
      hi <- if (!is.null(d$max)) d$max else beta_j
      function(x) punif(x, min = lo, max = hi)
    },

    "chisq" = {
      df <- d$df
      if (is.null(df)) stop("chisq distribution requires 'df'.")
      function(x) pchisq(x, df = df)
    },

    "t" = {
      df  <- d$df
      if (is.null(df)) stop("t distribution requires 'df'.")
      ncp <- d$ncp %||% 0
      function(x) pt(x, df = df, ncp = ncp)
    },

    "cauchy" = {
      loc <- d$location %||% 0;  sc <- d$scale %||% 1
      function(x) pcauchy(x, location = loc, scale = sc)
    },

    "logis" = {
      loc <- d$location %||% 0;  sc <- d$scale %||% 1
      function(x) plogis(x, location = loc, scale = sc)
    },

    "pareto" = {
      # Pareto: F(x) = 1 - (scale/x)^shape  for x >= scale
      sh <- d$shape;  sc <- d$scale
      if (is.null(sh) || is.null(sc))
        stop("pareto distribution requires 'shape' and 'scale'.")
      function(x) ifelse(x < sc, 0, 1 - (sc / x)^sh)
    },

    "triang" = {
      # Triangular(a, b, c) where c is the mode (c = b for right-triangular)
      a <- d$a;  b <- d$b;  cc <- d$c
      if (is.null(a) || is.null(b) || is.null(cc))
        stop("triang distribution requires 'a', 'b', and 'c' (mode).")
      .cdf_triang(a, b, cc)
    },

    stop("Unsupported distribution family: '", fam, "'. ",
         "Supported: norm, lnorm, gamma, exp, weibull, beta, unif, chisq, ",
         "t, cauchy, logis, pareto, triang.")
  )
}

#' Build a PDF function for a given distribution specification
#' Used by the analytic weight gradient (paper Section 4.2): dW_h/db involves
#' f_j(b) evaluated at the boundary.  Same family catalogue as .build_cdf.
#' @keywords internal
.build_pdf <- function(d, alpha_j, beta_j) {
  fam <- tolower(trimws(d$family))

  switch(fam,
    "norm"    = { mn <- d$mean %||% 0;  sd <- d$sd %||% 1
                  function(x) dnorm(x, mean = mn, sd = sd) },
    "lnorm"   = { ml <- d$meanlog %||% 0;  sl <- d$sdlog %||% 1
                  function(x) dlnorm(x, meanlog = ml, sdlog = sl) },
    "gamma"   = { sh <- d$shape
                  if (is.null(sh)) stop("gamma distribution requires 'shape'.")
                  rt <- if (!is.null(d$rate)) d$rate
                        else if (!is.null(d$scale)) 1 / d$scale else 1
                  function(x) dgamma(x, shape = sh, rate = rt) },
    "exp"     = { rt <- d$rate %||% 1
                  function(x) dexp(x, rate = rt) },
    "weibull" = { sh <- d$shape %||% 1;  sc <- d$scale %||% 1
                  function(x) dweibull(x, shape = sh, scale = sc) },
    "beta"    = { s1 <- d$shape1 %||% 1;  s2 <- d$shape2 %||% 1
                  rng <- beta_j - alpha_j
                  function(x) dbeta((x - alpha_j) / rng,
                                    shape1 = s1, shape2 = s2) / rng },
    "unif"    = { lo <- if (!is.null(d$min)) d$min else alpha_j
                  hi <- if (!is.null(d$max)) d$max else beta_j
                  function(x) dunif(x, min = lo, max = hi) },
    "chisq"   = { df <- d$df
                  if (is.null(df)) stop("chisq distribution requires 'df'.")
                  function(x) dchisq(x, df = df) },
    "t"       = { df <- d$df
                  if (is.null(df)) stop("t distribution requires 'df'.")
                  ncp <- d$ncp %||% 0
                  function(x) dt(x, df = df, ncp = ncp) },
    "cauchy"  = { loc <- d$location %||% 0;  sc <- d$scale %||% 1
                  function(x) dcauchy(x, location = loc, scale = sc) },
    "logis"   = { loc <- d$location %||% 0;  sc <- d$scale %||% 1
                  function(x) dlogis(x, location = loc, scale = sc) },
    "pareto"  = { sh <- d$shape;  sc <- d$scale
                  if (is.null(sh) || is.null(sc))
                    stop("pareto distribution requires 'shape' and 'scale'.")
                  function(x) ifelse(x < sc, 0, sh * sc^sh / x^(sh + 1)) },
    "triang"  = { a <- d$a;  b <- d$b;  cc <- d$c
                  if (is.null(a) || is.null(b) || is.null(cc))
                    stop("triang distribution requires 'a', 'b', and 'c' (mode).")
                  .pdf_triang(a, b, cc) },
    stop("Unsupported distribution family: '", fam, "'.")
  )
}

#' Build a quantile function for a given distribution specification
#' Used by the Algorithm 2 (MSAO) initialisation: random boundary starts are
#' drawn as \eqn{F_j^{-1}(u)} with \eqn{u} uniform, using the parametric quantile in
#' parametric mode.  Same family catalogue as .build_cdf.
#' @keywords internal
.build_quantile <- function(d, alpha_j, beta_j) {
  fam <- tolower(trimws(d$family))

  switch(fam,
    "norm"    = { mn <- d$mean %||% 0;  sd <- d$sd %||% 1
                  function(u) qnorm(u, mean = mn, sd = sd) },
    "lnorm"   = { ml <- d$meanlog %||% 0;  sl <- d$sdlog %||% 1
                  function(u) qlnorm(u, meanlog = ml, sdlog = sl) },
    "gamma"   = { sh <- d$shape
                  if (is.null(sh)) stop("gamma distribution requires 'shape'.")
                  rt <- if (!is.null(d$rate)) d$rate
                        else if (!is.null(d$scale)) 1 / d$scale else 1
                  function(u) qgamma(u, shape = sh, rate = rt) },
    "exp"     = { rt <- d$rate %||% 1
                  function(u) qexp(u, rate = rt) },
    "weibull" = { sh <- d$shape %||% 1;  sc <- d$scale %||% 1
                  function(u) qweibull(u, shape = sh, scale = sc) },
    "beta"    = { s1 <- d$shape1 %||% 1;  s2 <- d$shape2 %||% 1
                  rng <- beta_j - alpha_j
                  function(u) alpha_j + rng * qbeta(u, shape1 = s1, shape2 = s2) },
    "unif"    = { lo <- if (!is.null(d$min)) d$min else alpha_j
                  hi <- if (!is.null(d$max)) d$max else beta_j
                  function(u) qunif(u, min = lo, max = hi) },
    "chisq"   = { df <- d$df
                  if (is.null(df)) stop("chisq distribution requires 'df'.")
                  function(u) qchisq(u, df = df) },
    "t"       = { df <- d$df
                  if (is.null(df)) stop("t distribution requires 'df'.")
                  ncp <- d$ncp %||% 0
                  function(u) qt(u, df = df, ncp = ncp) },
    "cauchy"  = { loc <- d$location %||% 0;  sc <- d$scale %||% 1
                  function(u) qcauchy(u, location = loc, scale = sc) },
    "logis"   = { loc <- d$location %||% 0;  sc <- d$scale %||% 1
                  function(u) qlogis(u, location = loc, scale = sc) },
    "pareto"  = { sh <- d$shape;  sc <- d$scale
                  if (is.null(sh) || is.null(sc))
                    stop("pareto distribution requires 'shape' and 'scale'.")
                  function(u) sc * (1 - u)^(-1 / sh) },
    "triang"  = { a <- d$a;  b <- d$b;  cc <- d$c
                  if (is.null(a) || is.null(b) || is.null(cc))
                    stop("triang distribution requires 'a', 'b', and 'c' (mode).")
                  .q_triang(a, b, cc) },
    stop("Unsupported distribution family: '", fam, "'.")
  )
}

# Triangular PDF, vectorised
.pdf_triang <- function(a, b, cc) {
  function(x) {
    out <- numeric(length(x))
    mid_lo <- x >= a & x < cc
    mid_hi <- x >= cc & x <= b
    if (cc > a) out[mid_lo] <- 2 * (x[mid_lo] - a) / ((b - a) * (cc - a))
    if (b > cc) out[mid_hi] <- 2 * (b - x[mid_hi]) / ((b - a) * (b - cc))
    if (cc == b) out[x == b] <- 2 / (b - a)      # right-triangular apex
    out
  }
}

# Triangular quantile, vectorised (inverse of .cdf_triang)
.q_triang <- function(a, b, cc) {
  Fc <- if (b > a) (cc - a) / (b - a) else 0.5
  function(u) {
    u   <- pmin(pmax(u, 0), 1)
    out <- numeric(length(u))
    lo  <- u <= Fc
    out[lo]  <- a + sqrt(u[lo] * (b - a) * (cc - a))
    out[!lo] <- b - sqrt((1 - u[!lo]) * (b - a) * (b - cc))
    out
  }
}

# Triangular CDF, vectorised
.cdf_triang <- function(a, b, cc) {
  function(x) {
    out <- numeric(length(x))
    lo  <- x <= a
    hi  <- x >= b
    mid_lo <- !lo & x < cc
    mid_hi <- !hi & x >= cc & !lo
    out[lo]     <- 0
    out[hi]     <- 1
    out[mid_lo] <- (x[mid_lo] - a)^2 / ((b - a) * (cc - a))
    out[mid_hi] <- 1 - (b - x[mid_hi])^2 / ((b - a) * (b - cc))
    out
  }
}

# Minimal null-coalescing helper (avoids dependency on rlang)
`%||%` <- function(a, b) if (!is.null(a)) a else b
