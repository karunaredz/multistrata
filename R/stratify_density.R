# =============================================================================
# stratify_density.R
#
# Density-only (data-free) multivariate optimum stratification for the case in
# which the stratification variables ARE the study variables (Case 1 of
# Khan & Reddy).  No population data frame is required: the marginal densities
# alone determine every stratum quantity.
#
#   * Stratum weights   W_h    = prod_j [ F_j(b_hj) - F_j(b_{h-1,j}) ]   (as in
#     the parametric mode of stratify_multi()), computed from the marginal CDFs.
#   * Stratum means     mu_hk  and variances sigma2_hk of X_k are obtained by
#     numerical quadrature of the marginal density over each stratum interval:
#       mu   = ( \int x   f(x) dx ) / ( \int f(x) dx ),
#       s2   = ( \int x^2 f(x) dx ) / ( \int f(x) dx ) - mu^2,
#     evaluated with stats::integrate().  Because X_k is one of the
#     stratification variables and the margins are independent, mu_hk and
#     sigma2_hk depend only on the k-th coordinate interval of stratum h.
#
# The optimiser and objective are identical to stratify_multi(): a multi-start
# alternating scheme (derivative-free COBYLA boundary update + closed-form
# generalised Neyman allocation).  In the large-N limit the empirical
# (data-based) mode of stratify_multi() converges to this analytic mode; here
# there is no sampling noise and no dependence on a realised frame -- the exact
# analogue of the density-based LINGO formulation of Khan & Reddy.
#
# Reference:
#   Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#   Allocation for Multivariate Surveys (upcoming methods paper).
# =============================================================================


#' Data-Free Multivariate Optimum Stratification from Marginal Densities
#'
#' Determines optimum strata boundaries and sample allocation directly from the
#' marginal distributions of the stratification variables, with no population
#' data frame, for the case in which the stratification variables are
#' themselves the study variables (Case 1 of Khan & Reddy).  Stratum weights
#' are computed from the marginal CDFs and stratum means and variances by
#' quadrature of the marginal densities, so the design is exact up to
#' integration tolerance and free of sampling noise.  This is the analytic
#' analogue of \code{\link{stratify_multi}} and matches it in the large-\eqn{N}
#' limit.
#'
#' @param distributions A named list of length \eqn{p}: one marginal
#'   specification per stratification variable, in the same format as the
#'   \code{distributions} argument of \code{\link{stratify_multi}} (each a list
#'   with a \code{family} field and its parameters).  Supported families are
#'   those of \code{\link{compute_stratum_stats_parametric}}.
#' @param H_vec Integer vector of length \eqn{p}: number of strata per variable.
#' @param n_total Total sample size to allocate.
#' @param N Population size \eqn{N} (used only to form the finite-population
#'   corrections \eqn{N_h = W_h N}).
#' @param ranges Optional list of length \eqn{p}, each a length-2 numeric
#'   \code{c(lower, upper)} giving the support \eqn{[\alpha_j, \beta_j]} over
#'   which variable \eqn{j} is stratified.  When \code{NULL} (default) each
#'   range defaults to the family's \code{c(1e-6, 1 - 1e-6)} quantiles.  Supply
#'   explicit ranges to match a truncated support (e.g. \code{c(0, 10)} for an
#'   exponential capped at 10).
#' @param omega Study-variable weight vector of length \eqn{p} summing to one.
#'   Defaults to equal weights.  (Here \eqn{q = p}: study variable \eqn{k} is
#'   stratification variable \eqn{k}.)
#' @param constraint,cost_total,cost_fixed,cost_per_unit,cv_targets,n_min As in
#'   \code{\link{stratify_multi}} / \code{\link{neyman_alloc}}.
#' @param n_starts Number of random starts.  Default \code{20} (the analytic
#'   objective is smooth, so fewer starts are needed than in the empirical
#'   mode).
#' @param max_iter Maximum alternating-optimisation iterations per start.
#'   Default \code{100}.
#' @param tol_phi Relative objective convergence tolerance.  Default
#'   \code{1e-8}.
#' @param subdivisions Maximum subdivisions passed to \code{stats::integrate}
#'   for the stratum-moment quadrature.  Default \code{200}.
#' @param verbose If \code{TRUE}, prints per-start progress.
#'
#' @return An object of class \code{"multistrata"} (with
#'   \code{.meta$mode = "density"}), carrying \code{boundaries},
#'   \code{allocation}, \code{objective}, \code{optimality_gap},
#'   \code{stratum_stats}, \code{ranges}, \code{converged} and \code{call}.
#'   \code{print()} and \code{summary()} apply; \code{plot_strata()} does not
#'   (there is no data to scatter).
#'
#' @examples
#' \donttest{
#' # Khan & Reddy bivariate example: Exponential(1) x right-triangular(0,2),
#' # 4 x 3 = 12 strata, N = 2000, n = 700, no data at all.
#' res <- stratify_multi_density(
#'   distributions = list(x = list(family = "exp",    rate = 1),
#'                        y = list(family = "triang", a = 0, b = 2, c = 0)),
#'   H_vec   = c(4, 3),
#'   n_total = 700,
#'   N       = 2000,
#'   ranges  = list(c(0, 10), c(0, 2)))
#' print(res)
#' }
#'
#' @section Reproducibility:
#' As in \code{\link{stratify_multi}}, the multi-start scheme draws its initial
#' boundaries at random and sets no seed, so call
#' \code{\link[base]{set.seed}} beforehand when the design must be
#' reproducible.  The objective here is smooth rather than a step function, so
#' the spread across starts is typically smaller than in the empirical mode,
#' but it is not zero.
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{neyman_alloc}},
#'   \code{\link{compute_stratum_stats_parametric}}
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#' Allocation for Multivariate Surveys (upcoming methods paper).
#'
#' @export
stratify_multi_density <- function(distributions,
                                   H_vec,
                                   n_total,
                                   N,
                                   ranges        = NULL,
                                   omega         = NULL,
                                   constraint    = c("n_total", "cost"),
                                   cost_total    = NULL,
                                   cost_fixed    = 0,
                                   cost_per_unit = 1,
                                   cv_targets    = NULL,
                                   n_min         = 2L,
                                   n_starts      = 20L,
                                   max_iter      = 100L,
                                   tol_phi       = 1e-8,
                                   subdivisions  = 200L,
                                   verbose       = FALSE) {

  mc         <- match.call()
  constraint <- match.arg(constraint)

  # ------------------------------------------------------------------
  # 1.  Validate inputs
  # ------------------------------------------------------------------
  if (!is.list(distributions))
    stop("'distributions' must be a list of marginal specifications.")
  p <- length(distributions)
  if (length(H_vec) != p)
    stop("'H_vec' must have one entry per distribution (length p = ", p, ").")
  if (any(H_vec < 2L))
    stop("Each element of 'H_vec' must be >= 2.")
  for (j in seq_len(p))
    if (!is.list(distributions[[j]]) || is.null(distributions[[j]]$family))
      stop("distributions[[", j, "]] must be a list with a 'family' field.")

  q <- p                                  # study variables = stratification vars
  if (is.null(omega)) omega <- rep(1 / q, q)
  if (length(omega) != q) stop("'omega' must have length p = ", p, ".")
  if (abs(sum(omega) - 1) > 1e-8) stop("'omega' must sum to 1.")
  if (any(omega <= 0)) stop("All elements of 'omega' must be strictly positive.")
  if (!is.null(cv_targets) && length(cv_targets) != q)
    stop("'cv_targets' must have length p = ", p, ".")

  x_names <- if (!is.null(names(distributions))) names(distributions)
             else paste0("X", seq_len(p))

  # Resolve ranges (support per variable)
  if (is.null(ranges)) {
    ranges <- lapply(seq_len(p), function(j) {
      qf <- .build_quantile(distributions[[j]], 0, 1)
      c(qf(1e-6), qf(1 - 1e-6))
    })
  }
  if (length(ranges) != p)
    stop("'ranges' must be a list of length p = ", p, ".")
  alpha <- vapply(ranges, function(r) r[1L], numeric(1))
  beta  <- vapply(ranges, function(r) r[2L], numeric(1))
  if (any(!is.finite(alpha)) || any(!is.finite(beta)) || any(beta <= alpha))
    stop("Each range must be finite with lower < upper. ",
         "Supply 'ranges' explicitly for heavy-tailed families.")

  H_total <- prod(H_vec)
  qfuns   <- lapply(seq_len(p), function(j)
                    .build_quantile(distributions[[j]], alpha[j], beta[j]))

  # ------------------------------------------------------------------
  # 2.  Multi-start alternating optimisation
  # ------------------------------------------------------------------
  dummy    <- matrix(0, nrow = 1L, ncol = p)     # only ncol() is used
  best_phi <- Inf; best_B <- NULL; best_alloc <- NULL; best_conv <- FALSE

  for (r in seq_len(n_starts)) {
    if (verbose) cat(sprintf("  Start %d / %d ...\n", r, n_starts))
    B_init <- .random_init_boundaries(dummy, H_vec, alpha, beta, qfuns)
    marg0  <- .density_margins(distributions, B_init, H_vec, ranges, subdivisions)
    asm0   <- .assemble_from_margins(marg0, H_vec, N)
    n_init <- neyman_alloc(sigma2 = asm0$sigma2, W_h = asm0$W_h, omega = omega,
                           n_total = n_total, N_h = asm0$N_h, n_min = n_min,
                           constraint = constraint, cost_total = cost_total,
                           cost_fixed = cost_fixed, cost_per_unit = cost_per_unit,
                           cv_targets = cv_targets, mu = asm0$mu)

    ao <- .ao_core_density(
      distributions, B_init, n_init, H_vec, ranges, alpha, beta, omega, N,
      n_total, constraint, cost_total, cost_fixed, cost_per_unit,
      cv_targets, n_min, max_iter, tol_phi, subdivisions)

    if (ao$phi < best_phi) {
      best_phi <- ao$phi; best_B <- ao$B
      best_alloc <- ao$n_alloc; best_conv <- ao$converged
      if (verbose) cat(sprintf("    -> New best Phi = %.8g\n", best_phi))
    }
  }

  # ------------------------------------------------------------------
  # 3.  Post-processing
  # ------------------------------------------------------------------
  marg <- .density_margins(distributions, best_B, H_vec, ranges, subdivisions)
  asm  <- .assemble_from_margins(marg, H_vec, N)
  labels <- .stratum_index_labels(H_vec)

  ss <- data.frame(stratum = labels, N_h = asm$N_h, W_h = asm$W_h,
                   n_h = as.integer(round(best_alloc)),
                   stringsAsFactors = FALSE)
  for (k in seq_len(q)) {
    ss[[paste0("mu_k",     k)]] <- asm$mu[, k]
    ss[[paste0("sigma2_k", k)]] <- asm$sigma2[, k]
  }
  ss[["A_h"]] <- .compute_A_h(asm$sigma2, asm$W_h, omega)

  gap <- {
    lb <- variance_lower_bound(asm$sigma2, asm$W_h, asm$N_h, omega, n_total)
    if (lb < 1e-15) 0 else (best_phi - lb) / abs(lb)
  }
  cv_ach <- attr(best_alloc, "cv_achieved")

  structure(
    list(
      boundaries     = best_B,
      allocation     = best_alloc,
      objective      = best_phi,
      optimality_gap = gap,
      kkt_residual   = NULL,
      kkt_ok         = NA,
      stratum_stats  = ss,
      cv_achieved    = cv_ach,
      ranges         = ranges,
      n_starts_run   = n_starts,
      converged      = best_conv,
      call           = mc,
      .meta = list(
        p = p, q = q, H_vec = H_vec, H_total = H_total,
        N = N, n_total = n_total, omega = omega, constraint = constraint,
        cv_targets = cv_targets, method = "cobyla",
        alpha = alpha, beta = beta,
        x_names = x_names, y_names = x_names,
        mode = "density",
        distributions = distributions,
        auto_fit = NULL
      )
    ),
    class = "multistrata"
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# One marginal's per-interval probability, mean and variance by quadrature.
.density_one_margin <- function(d, cuts, range, Hj, subdivisions = 200L) {
  lo <- range[1L]; hi <- range[2L]
  cdf <- .build_cdf(d, lo, hi)
  pdf <- .build_pdf(d, lo, hi)
  full <- c(lo, cuts, hi)

  pj <- numeric(Hj); mj <- numeric(Hj); vj <- numeric(Hj)
  for (i in seq_len(Hj)) {
    a <- full[i]; b <- full[i + 1L]
    pj[i] <- max(cdf(b) - cdf(a), 0)
    if (b <= a || pj[i] <= 1e-12) {
      mj[i] <- (a + b) / 2; vj[i] <- 0; next
    }
    m1 <- tryCatch(stats::integrate(function(x) x * pdf(x), a, b,
                     subdivisions = subdivisions, rel.tol = 1e-8,
                     stop.on.error = FALSE)$value, error = function(e) NA_real_)
    m2 <- tryCatch(stats::integrate(function(x) x * x * pdf(x), a, b,
                     subdivisions = subdivisions, rel.tol = 1e-8,
                     stop.on.error = FALSE)$value, error = function(e) NA_real_)
    if (!is.finite(m1) || !is.finite(m2)) {
      mj[i] <- (a + b) / 2; vj[i] <- (b - a)^2 / 12; next
    }
    mean_i <- m1 / pj[i]
    mj[i]  <- mean_i
    vj[i]  <- max(m2 / pj[i] - mean_i^2, 0)
  }
  s <- sum(pj); if (s < 1e-12) s <- 1
  list(p = pj / s, mean = mj, var = vj)
}

# All marginals.
.density_margins <- function(distributions, B, H_vec, ranges, subdivisions = 200L) {
  lapply(seq_along(H_vec), function(j)
    .density_one_margin(distributions[[j]], B[[j]], ranges[[j]],
                        H_vec[j], subdivisions))
}

# Assemble stratum-level W_h, N_h, mu, sigma2 from the per-margin statistics.
.assemble_from_margins <- function(margins, H_vec, N) {
  p       <- length(H_vec)
  H_total <- prod(H_vec)
  grids   <- lapply(seq_len(p), function(j) seq_len(H_vec[j]))
  combos  <- as.matrix(expand.grid(grids, KEEP.OUT.ATTRS = FALSE))

  W_h <- apply(combos, 1L, function(rw)
    prod(vapply(seq_len(p), function(j) margins[[j]]$p[rw[j]], numeric(1))))
  W_h <- W_h / sum(W_h)
  N_h <- pmax(round(W_h * N), 1L)

  mu     <- matrix(0, H_total, p)
  sigma2 <- matrix(0, H_total, p)
  for (k in seq_len(p)) {
    mu[, k]     <- margins[[k]]$mean[combos[, k]]
    sigma2[, k] <- margins[[k]]$var[ combos[, k]]
  }
  list(W_h = W_h, N_h = N_h, mu = mu, sigma2 = sigma2)
}

# One alternating-optimisation run (density objective).
.ao_core_density <- function(distributions, B, n_alloc, H_vec, ranges,
                             alpha, beta, omega, N, n_total, constraint,
                             cost_total, cost_fixed, cost_per_unit, cv_targets,
                             n_min, max_iter, tol_phi, subdivisions) {
  p <- length(H_vec)

  marg <- .density_margins(distributions, B, H_vec, ranges, subdivisions)
  asm  <- .assemble_from_margins(marg, H_vec, N)
  phi_prev  <- .composite_obj(asm$sigma2, asm$W_h, n_alloc, asm$N_h, omega)
  phi_hist  <- phi_prev
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # (A) Boundary update, one block per variable (COBYLA)
    for (j in seq_len(p)) {
      b_j_init <- B[[j]]; nb <- length(b_j_init)
      lbj   <- rep(alpha[j] + 1e-10, nb)
      ubj   <- rep(beta[j]  - 1e-10, nb)
      delta <- 1e-6 * (beta[j] - alpha[j])

      obj_fun <- function(b_j) {
        mj <- .density_one_margin(distributions[[j]], sort(b_j), ranges[[j]],
                                  H_vec[j], subdivisions)
        mt <- marg; mt[[j]] <- mj
        a  <- .assemble_from_margins(mt, H_vec, N)
        .composite_obj(a$sigma2, a$W_h, n_alloc, a$N_h, omega)
      }
      ineq_fun <- if (nb > 1L) function(b_j) delta - diff(b_j) else NULL

      phi_before <- .composite_obj(asm$sigma2, asm$W_h, n_alloc, asm$N_h, omega)
      res <- tryCatch(
        nloptr::nloptr(x0 = b_j_init, eval_f = obj_fun, lb = lbj, ub = ubj,
                       eval_g_ineq = ineq_fun,
                       opts = list(algorithm = "NLOPT_LN_COBYLA",
                                   xtol_rel = 1e-6, ftol_rel = 1e-8,
                                   maxeval = 300L)),
        error = function(e) NULL)

      if (!is.null(res) && res$objective < phi_before - 1e-14) {
        B[[j]]   <- sort(res$solution)
        marg[[j]] <- .density_one_margin(distributions[[j]], B[[j]], ranges[[j]],
                                         H_vec[j], subdivisions)
        asm      <- .assemble_from_margins(marg, H_vec, N)
      }
    }

    # (B) Allocation update
    n_alloc <- neyman_alloc(
      sigma2 = asm$sigma2, W_h = asm$W_h, omega = omega, n_total = n_total,
      N_h = asm$N_h, n_min = n_min, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit, cv_targets = cv_targets, mu = asm$mu)

    phi_new  <- .composite_obj(asm$sigma2, asm$W_h, n_alloc, asm$N_h, omega)
    phi_hist <- c(phi_hist, phi_new)
    if (abs(phi_new - phi_prev) / (1 + abs(phi_prev)) < tol_phi) {
      converged <- TRUE; phi_prev <- phi_new; break
    }
    phi_prev <- phi_new
  }

  list(B = B, n_alloc = n_alloc, phi = phi_prev, converged = converged,
       phi_trace = phi_hist)
}
