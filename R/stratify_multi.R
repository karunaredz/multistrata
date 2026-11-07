# =============================================================================
# stratify_multi.R
#
# Core function: simultaneous multivariate optimum stratification and sample
# allocation via multi-start alternating optimisation (MSAO).
#
# Fixes and additions in this revision (relative to the 0.2.0 sources):
#   1. .compute_stratum_stats(): A_h now computed with the corrected
#      .compute_A_h() (see neyman_alloc.R Fix 1) -- the old sweep() call
#      misaligned omega/W_h whenever H != q.
#   2. COBYLA ordering constraint direction corrected.  nloptr requires
#      g(x) <= 0; the old code returned diff(b) - delta, which *enforces*
#      diff(b) <= delta (near-coincident boundaries) instead of
#      diff(b) >= delta.  Now returns delta - diff(b).
#   3. Multi-start initialisation now follows Algorithm 2 (MSAO) of the
#      paper: uniform order-statistic levels mapped through the quantile
#      function -- the parametric quantile when 'distributions' is given,
#      the empirical quantile otherwise -- plus jitter.  This makes the
#      random starts dense in the boundary space, as required by the
#      asymptotic global-optimality result (Proposition 2).
#   4. 'method' argument: boundary sub-problem solved by derivative-free
#      COBYLA (default) or by gradient-based SLSQP ("slsqp", parametric
#      mode only) using the closed-form weight gradient dW_h/db of
#      Section 4.2 of the paper; the empirical Y-moments are piecewise constant in B,
#      so their a.e.-zero derivative is omitted from the gradient.
#   5. tol_kkt is now used: a scaled KKT residual (allocation-stationarity
#      spread + boundary directional derivatives) is computed at the final
#      solution and reported; 'kkt_ok' flags whether it is below tol_kkt.
#   6. Cost-mode and CV-constrained allocation are passed through to the
#      corrected neyman_alloc() ('cost_fixed', 'cv_targets').
#   7. .optimality_gap() now delegates to variance_lower_bound() so the
#      reported gap and the exported bound cannot drift apart.
#
# Reference:
#   Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#   Allocation for Multivariate Surveys (upcoming methods paper).
# =============================================================================


#' Simultaneous Multivariate Optimum Stratification and Sample Allocation
#'
#' Jointly determines strata boundaries and sample allocation across multiple
#' stratification variables and multiple study variables by minimising a
#' weighted composite variance objective.  The algorithm alternates between a
#' boundary update (derivative-free COBYLA, or gradient-based SLSQP under the
#' parametric mode) and a closed-form generalised Neyman allocation update
#' (Theorem 1 in Reddy & Khan, forthcoming), embedded within a multi-start wrapper
#' to mitigate convergence to poor local minima.
#'
#' @param data A \code{data.frame} or \code{matrix} whose columns are the
#'   stratification (auxiliary) variables; rows represent population units.
#' @param study_vars A \code{data.frame} or \code{matrix} of study (target)
#'   variables with the same number of rows as \code{data}.
#' @param H_vec Integer vector of length \eqn{p} giving the number of strata
#'   for each stratification variable.  The total number of strata is
#'   \code{prod(H_vec)}.
#' @param n_total Total sample size to be allocated across all strata.
#' @param omega Weight vector of length \eqn{q} summing to one.  Defaults to
#'   equal weights.
#' @param constraint One of \code{"n_total"} (default) or \code{"cost"}.
#' @param cost_total Scalar budget \eqn{C^*} (used when
#'   \code{constraint = "cost"}).
#' @param cost_fixed Fixed overhead \eqn{c_0} subtracted from the budget.
#'   Default \code{0}.
#' @param cost_per_unit Scalar or length-\eqn{H} numeric vector of per-unit
#'   costs.  Default \code{1}.
#' @param cv_targets Optional numeric vector of length \eqn{q}: maximum
#'   coefficients of variation \eqn{\gamma_k} for the study-variable mean
#'   estimators (constraint (C3) in the paper).  When supplied, the
#'   allocation step solves the convex CV-constrained sub-problem instead of
#'   the closed form.  Default \code{NULL}.
#' @param n_min Minimum sample size per stratum.  Default \code{2}.
#' @param n_starts Number of random starts.  Default \code{50}.  See the
#'   \strong{Reproducibility} section: the starts are drawn at random, so set a
#'   seed when the design must be reproducible.
#' @param max_iter Maximum AO iterations per start.  Default \code{200}.
#' @param tol_phi Relative objective convergence tolerance.  Default
#'   \code{1e-8}.
#' @param tol_kkt Threshold for the scaled KKT residual reported at the
#'   final solution (secondary diagnostic).  Default \code{1e-6}.
#' @param method Boundary sub-problem solver: \code{"cobyla"} (default,
#'   derivative-free) or \code{"slsqp"} (gradient-based, requires
#'   \code{distributions}; uses the closed-form weight gradient of Section 4.2 of the paper).
#' @param verbose If \code{TRUE}, prints per-start progress.
#' @param trace If \code{TRUE}, the returned object carries \code{phi_trace},
#'   the per-iteration composite-variance sequence of the best start (used, for
#'   example, to visualise the monotone descent of Section 4.2).  Default
#'   \code{FALSE}.
#' @param distributions Controls the parametric mode.  Either \code{NULL}
#'   (default, empirical mode), the string \code{"auto"} (automatic
#'   AIC-based family selection per variable via
#'   \code{\link{fit_distributions}}, using the FitVerse package as an
#'   additional fitting engine when installed), or a named list of length
#'   \eqn{p} in which each element is a list with a \code{family} field and
#'   distribution parameters matching the base-R parameterisation (e.g.,
#'   \code{list(family = "lnorm", meanlog = 1.2, sdlog = 0.8)}).  Supported
#'   families: \code{"norm"}, \code{"lnorm"}, \code{"gamma"}, \code{"exp"},
#'   \code{"weibull"}, \code{"beta"}, \code{"unif"}, \code{"chisq"},
#'   \code{"t"}, \code{"cauchy"}, \code{"logis"}, \code{"pareto"},
#'   \code{"triang"}.  When supplied (or auto-selected), stratum weights
#'   \eqn{W_h} are computed analytically from marginal CDF differences
#'   rather than empirically from sample proportions.  Study-variable
#'   moments (\eqn{\sigma^2_{hk}}) are always estimated empirically.
#'
#' @return An object of class \code{"multistrata"} with components
#'   \code{boundaries}, \code{allocation}, \code{objective},
#'   \code{optimality_gap}, \code{kkt_residual} (the larger of the two blocks),
#'   \code{kkt_alloc} and \code{kkt_boundary} (the allocation and boundary KKT
#'   residuals reported separately), \code{kkt_ok}, \code{stratum_stats},
#'   \code{cv_achieved} (when \code{cv_targets} was supplied),
#'   \code{n_starts_run}, \code{n_starts_failed}, \code{objective_by_start} and
#'   \code{converged_by_start} (the objective reached and convergence flag for
#'   each start), \code{converged}, \code{phi_trace} (when \code{trace = TRUE}),
#'   and \code{call}.
#'
#' @examples
#' set.seed(42)
#' pop <- sim_population(N = 2000, p = 2, q = 2, scenario = "independent")
#' result <- stratify_multi(
#'   data       = pop[, c("X1", "X2")],
#'   study_vars = pop[, c("Y1", "Y2")],
#'   H_vec      = c(3, 3),
#'   n_total    = 500,
#'   omega      = c(0.5, 0.5),
#'   n_starts   = 5,
#'   verbose    = FALSE
#' )
#' print(result)
#' summary(result)
#'
#' @section Reproducibility:
#' The multi-start scheme draws its initial boundaries at random and does not
#' set a seed, so repeated calls on the same data can return slightly
#' different designs.  Call \code{\link[base]{set.seed}} beforehand whenever
#' the design has to be reproducible --- which includes any design you intend
#' to report:
#'
#' \preformatted{  set.seed(1)
#'   res <- stratify_multi(X, Y, H_vec = c(4, 3), n_total = 600)}
#'
#' The spread across starts narrows quickly as \code{n_starts} grows but does
#' not vanish.  On a 12,000-unit bivariate population it was roughly 14\% in
#' \eqn{\Phi^*} at \code{n_starts = 1} and about 0.1\% at
#' \code{n_starts = 50}.  Treat a large spread as a signal to raise
#' \code{n_starts} rather than as noise to be re-rolled, and compare
#' \code{optimality_gap} across runs when in doubt.
#'
#' @seealso \code{\link{neyman_alloc}}, \code{\link{plot_strata}},
#'   \code{\link{benchmark_designs}}
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#' Allocation for Multivariate Surveys (upcoming methods paper).
#'
#' @export
stratify_multi <- function(data,
                           study_vars,
                           H_vec,
                           n_total,
                           omega         = NULL,
                           constraint    = c("n_total", "cost"),
                           cost_total    = NULL,
                           cost_fixed    = 0,
                           cost_per_unit = 1,
                           cv_targets    = NULL,
                           n_min         = 2L,
                           n_starts      = 50L,
                           max_iter      = 200L,
                           tol_phi       = 1e-8,
                           tol_kkt       = 1e-6,
                           method        = c("cobyla", "slsqp"),
                           verbose       = FALSE,
                           distributions = NULL,
                           trace         = FALSE) {

  mc         <- match.call()
  constraint <- match.arg(constraint)
  method     <- match.arg(method)

  # ------------------------------------------------------------------
  # 1.  Input validation
  # ------------------------------------------------------------------
  data       <- as.matrix(data)
  study_vars <- as.matrix(study_vars)

  N <- nrow(data);  p <- ncol(data);  q <- ncol(study_vars)

  if (nrow(study_vars) != N)
    stop("'data' and 'study_vars' must have the same number of rows.")
  if (length(H_vec) != p)
    stop("'H_vec' must have one entry per column of 'data'.")
  if (any(H_vec < 2L))
    stop("Each element of 'H_vec' must be >= 2.")

  if (is.null(omega)) omega <- rep(1 / q, q)
  if (length(omega) != q)
    stop("'omega' must have length equal to ncol(study_vars).")
  if (abs(sum(omega) - 1) > 1e-8)
    stop("'omega' must sum to 1.")
  if (any(omega <= 0))
    stop("All elements of 'omega' must be strictly positive.")
  if (!is.null(cv_targets) && length(cv_targets) != q)
    stop("'cv_targets' must have length equal to ncol(study_vars).")

  # Auto-parametric mode: select marginal families by AIC (FitVerse-aware)
  auto_fit <- NULL
  if (identical(distributions, "auto")) {
    if (verbose) message("  [auto] selecting marginal families by AIC ...")
    distributions <- fit_distributions(data, verbose = verbose)
    auto_fit <- attr(distributions, "fit_table")
    attr(distributions, "fit_table") <- NULL
    attr(distributions, "chosen")    <- NULL
  }

  # Validate distributions list if supplied
  parametric <- !is.null(distributions)
  if (parametric) {
    if (!is.list(distributions) || length(distributions) != p)
      stop("'distributions' must be a list of length p = ", p,
           " (one entry per stratification variable).")
    for (j in seq_len(p)) {
      if (!is.list(distributions[[j]]) || is.null(distributions[[j]]$family))
        stop("distributions[[", j, "]] must be a list with at least a 'family' field.")
    }
    if (verbose)
      cat("  [parametric mode] using analytical stratum weights from",
          paste(sapply(distributions, `[[`, "family"), collapse = ", "), "\n")
  }
  if (method == "slsqp" && !parametric)
    stop("method = 'slsqp' requires 'distributions' (parametric mode): ",
         "the analytic weight gradient of eq. (17) is only available when ",
         "the marginal CDFs are specified.  Use method = 'cobyla' for the ",
         "empirical mode.")

  H_total <- prod(H_vec)
  if (n_total < H_total * n_min)
    warning(sprintf(
      "n_total (%d) < H_total x n_min (%d x %d); some strata may be under-allocated.",
      n_total, H_total, n_min))

  # Support range + names, kept in .meta for methods
  alpha   <- apply(data, 2L, min)
  beta    <- apply(data, 2L, max)
  x_names <- if (!is.null(colnames(data)))       colnames(data)       else paste0("X", seq_len(p))
  y_names <- if (!is.null(colnames(study_vars))) colnames(study_vars) else paste0("Y", seq_len(q))

  # Pre-build quantile functions for Algorithm 2 initialisation
  qfuns <- if (parametric) {
    lapply(seq_len(p), function(j)
      .build_quantile(distributions[[j]], alpha[j], beta[j]))
  } else NULL

  # ------------------------------------------------------------------
  # 2.  Multi-start alternating optimisation (Algorithm 2, MSAO)
  # ------------------------------------------------------------------
  best_phi <- Inf;  best_B <- NULL;  best_alloc <- NULL;  best_conv <- FALSE
  best_trace <- NULL
  phi_by_start  <- rep(NA_real_, n_starts)   # objective reached by each start
  conv_by_start <- rep(NA, n_starts)         # did each start converge?
  n_fail        <- 0L                        # starts that failed to run

  for (r in seq_len(n_starts)) {
    if (verbose) cat(sprintf("  Start %d / %d ...\n", r, n_starts))

    B_init <- .random_init_boundaries(data, H_vec, alpha, beta, qfuns)
    n_init <- .proportional_alloc(data, B_init, H_vec, n_total, n_min)

    ao_res <- tryCatch(
      .ao_core(
        data = data, study_vars = study_vars, B = B_init, n_alloc = n_init,
        H_vec = H_vec, alpha = alpha, beta = beta, omega = omega,
        n_total = n_total, constraint = constraint, cost_total = cost_total,
        cost_fixed = cost_fixed, cost_per_unit = cost_per_unit,
        cv_targets = cv_targets, n_min = n_min,
        max_iter = max_iter, tol_phi = tol_phi, method = method,
        distributions = distributions
      ),
      error = function(e) NULL
    )

    if (is.null(ao_res)) { n_fail <- n_fail + 1L; next }

    phi_by_start[r]  <- ao_res$phi
    conv_by_start[r] <- ao_res$converged

    if (ao_res$phi < best_phi) {
      best_phi <- ao_res$phi;  best_B <- ao_res$B
      best_alloc <- ao_res$n_alloc;  best_conv <- ao_res$converged
      best_trace <- ao_res$phi_trace
      if (verbose) cat(sprintf("    -> New best Phi = %.8f\n", best_phi))
    }
  }

  if (is.null(best_B))
    stop(sprintf(
      "All %d optimisation start(s) failed; check the inputs, distributions, or constraints.",
      n_starts))

  # ------------------------------------------------------------------
  # 3.  Post-processing
  # ------------------------------------------------------------------
  strat_stats <- .compute_stratum_stats(
    data, study_vars, best_B, best_alloc, H_vec, omega, N, x_names, y_names,
    distributions = distributions, alpha = alpha, beta = beta)
  gap <- .optimality_gap(strat_stats, best_phi, n_total, omega)

  kkt <- .kkt_residual(
    data, study_vars, best_B, best_alloc, H_vec, alpha, beta, omega,
    distributions, n_min)

  cv_ach <- attr(best_alloc, "cv_achieved")
  if (!is.null(cv_targets) && !isTRUE(attr(best_alloc, "cv_feasible")))
    warning(sprintf(
      "CV targets were not attained by the final design (achieved: %s).",
      paste(sprintf("%.4f", cv_ach), collapse = ", ")))

  structure(
    list(
      boundaries     = best_B,
      allocation     = best_alloc,
      objective      = best_phi,
      optimality_gap = gap,
      kkt_residual   = kkt$max,
      kkt_alloc      = kkt$alloc,
      kkt_boundary   = kkt$boundary,
      kkt_ok         = kkt$max < tol_kkt,
      stratum_stats  = strat_stats,
      cv_achieved    = cv_ach,
      n_starts_run   = n_starts,
      n_starts_failed = n_fail,
      objective_by_start = phi_by_start,
      converged_by_start = conv_by_start,
      converged      = best_conv,
      phi_trace      = if (isTRUE(trace)) best_trace else NULL,
      call           = mc,
      .meta = list(
        p = p, q = q, H_vec = H_vec, H_total = H_total,
        N = N, n_total = n_total, omega = omega, constraint = constraint,
        cv_targets = cv_targets, method = method,
        alpha = alpha, beta = beta,
        x_names = x_names, y_names = y_names,
        mode = if (parametric) "parametric" else "empirical",
        distributions = distributions,
        auto_fit = auto_fit
      )
    ),
    class = "multistrata"
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# Algorithm 2, lines 3-5: draw H_j - 1 uniform order-statistic levels, map
# them through the quantile function (parametric F_j^{-1} when available,
# empirical otherwise), add jitter, clip and enforce strict ordering.
# FIX 3: the old code used *fixed* equally spaced quantiles + jitter, which
# is not dense in the boundary space and therefore not covered by the
# asymptotic global-optimality argument (Proposition 2 of the paper).
.random_init_boundaries <- function(data, H_vec, alpha, beta, qfuns = NULL) {
  p <- ncol(data)
  B <- vector("list", p)
  for (j in seq_len(p)) {
    H_j <- H_vec[j]
    nb  <- H_j - 1L
    u   <- sort(runif(nb, min = 0.02, max = 0.98))
    qpts <- if (!is.null(qfuns)) qfuns[[j]](u)
            else quantile(data[, j], probs = u, names = FALSE)
    zeta   <- (beta[j] - alpha[j]) / (4 * H_j)
    jitter <- runif(nb, min = -zeta, max = zeta)
    bpts   <- sort(pmax(alpha[j] + 1e-10,
                        pmin(beta[j] - 1e-10, qpts + jitter)))
    # Enforce strictly increasing, well-separated points
    min_gap <- 1e-6 * (beta[j] - alpha[j])
    if (nb > 1L) {
      for (i in 2:nb)
        if (bpts[i] - bpts[i - 1L] < min_gap)
          bpts[i] <- bpts[i - 1L] + min_gap
      bpts <- pmin(bpts, beta[j] - 1e-10)
    }
    if (any(!is.finite(bpts)) || (nb > 1L && any(diff(bpts) <= 0)))
      bpts <- alpha[j] + (beta[j] - alpha[j]) * seq_len(nb) / H_j
    B[[j]] <- bpts
  }
  B
}

.proportional_alloc <- function(data, B, H_vec, n_total, n_min) {
  N          <- nrow(data)
  H_total    <- prod(H_vec)
  strat_id   <- .assign_strata(data, B, H_vec)
  counts     <- tabulate(strat_id, nbins = H_total)
  raw        <- pmax(n_min, counts / N * n_total)
  alloc      <- .round_to_total(raw, as.integer(round(n_total)), n_min,
                                pmax(counts, n_min))
  setNames(as.numeric(alloc), .stratum_index_labels(H_vec))
}

.assign_strata <- function(data, B, H_vec) {
  N          <- nrow(data);  p <- ncol(data)
  idx_matrix <- matrix(0L, nrow = N, ncol = p)
  for (j in seq_len(p)) {
    breaks          <- c(-Inf, B[[j]], Inf)
    idx_j           <- as.integer(cut(data[, j], breaks = breaks,
                                      right = FALSE, labels = FALSE)) - 1L
    idx_matrix[, j] <- pmax(0L, pmin(H_vec[j] - 1L, idx_j))
  }
  # Mixed-radix encoding: stratum h = 1 + i1*1 + i2*H1 + i3*H1*H2 + ...
  multiplier <- cumprod(c(1L, H_vec[-length(H_vec)]))
  1L + as.integer(idx_matrix %*% multiplier)
}

.stratum_index_labels <- function(H_vec) {
  p      <- length(H_vec)
  grids  <- lapply(seq_len(p), function(j) seq_len(H_vec[j]))
  combos <- expand.grid(grids, KEEP.OUT.ATTRS = FALSE)
  apply(combos, 1L, function(row) paste0("h", seq_len(p), ".", row, collapse = "_"))
}

# Empirical stratum means/variances (study vars) + empirical or parametric W_h.
.compute_stratum_means_vars <- function(data, study_vars, B, H_vec,
                                        distributions = NULL,
                                        alpha = NULL, beta = NULL) {
  N        <- nrow(data);  q <- ncol(study_vars);  H_total <- prod(H_vec)
  strat_id <- .assign_strata(data, B, H_vec)
  counts   <- tabulate(strat_id, nbins = H_total)

  mu     <- matrix(0, nrow = H_total, ncol = q)
  sigma2 <- matrix(0, nrow = H_total, ncol = q)
  for (h in seq_len(H_total)) {
    idx <- which(strat_id == h)
    if (length(idx) < 2L) {
      mu[h, ]     <- colMeans(study_vars)
      sigma2[h, ] <- apply(study_vars, 2L, var)
    } else {
      ys          <- study_vars[idx, , drop = FALSE]
      mu[h, ]     <- colMeans(ys)
      sigma2[h, ] <- apply(ys, 2L, var)
    }
  }

  if (!is.null(distributions)) {
    ps <- compute_stratum_stats_parametric(
      data       = data,
      study_vars = study_vars,
      B          = B,
      H_vec      = H_vec,
      alpha      = if (!is.null(alpha)) alpha else apply(data, 2L, min),
      beta       = if (!is.null(beta))  beta  else apply(data, 2L, max),
      dist_list  = distributions
    )
    W_h_use <- ps$W_h
    N_h_use <- ps$N_h
  } else {
    W_h_use <- counts / N
    N_h_use <- counts
  }

  list(mu       = mu,
       sigma2   = sigma2,
       counts   = counts,
       strat_id = strat_id,
       W_h      = W_h_use,
       N_h      = N_h_use)
}

.composite_obj <- function(sigma2_mat, W_h, n_alloc, N_h, omega) {
  fpc <- pmax(1 / pmax(n_alloc, 1) - 1 / pmax(N_h, 1), 0)
  sum(omega * colSums(W_h^2 * sigma2_mat * fpc))
}


# .ao_core --------------------------------------------------------------------
# One AO run (Algorithm 1): boundary update (COBYLA or SLSQP) alternating
# with the generalised Neyman / CV-constrained allocation update.
.ao_core <- function(data, study_vars, B, n_alloc, H_vec, alpha, beta,
                     omega, n_total, constraint, cost_total, cost_fixed,
                     cost_per_unit, cv_targets, n_min, max_iter, tol_phi,
                     method, distributions = NULL) {

  N <- nrow(data);  p <- ncol(data)

  sv        <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                           distributions, alpha, beta)
  W_h       <- sv$W_h
  N_h       <- sv$N_h
  phi_prev  <- .composite_obj(sv$sigma2, W_h, n_alloc, N_h, omega)
  phi_hist  <- phi_prev            # per-iteration objective trace
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # (A) Boundary update, one block per stratification variable
    for (j in seq_len(p)) {
      b_j_init     <- B[[j]];  nb <- length(b_j_init)
      lb           <- rep(alpha[j] + 1e-10, nb)
      ub           <- rep(beta[j]  - 1e-10, nb)
      delta        <- 1e-6 * (beta[j] - alpha[j])
      phi_before_j <- .composite_obj(sv$sigma2, W_h, n_alloc, N_h, omega)

      obj_fun <- function(b_j) {
        B_try  <- B;  B_try[[j]] <- sort(b_j)
        sv_try <- .compute_stratum_means_vars(data, study_vars, B_try, H_vec,
                                              distributions, alpha, beta)
        .composite_obj(sv_try$sigma2, sv_try$W_h, n_alloc, sv_try$N_h, omega)
      }

      # FIX 2: nloptr convention is g(x) <= 0.  Minimum-spacing constraint
      # b_{i+1} - b_i >= delta  must be written  delta - diff(b) <= 0.
      # (The old code returned diff(b) - delta, i.e. diff(b) <= delta,
      # which instructs the solver to *collapse* the boundaries.)
      ineq_fun <- if (nb > 1L) function(b_j) delta - diff(b_j) else NULL

      res_nlopt <-
        if (method == "slsqp") {
          grad_fun <- function(b_j)
            .phi_grad_boundary(b_j, j, B, data, study_vars, H_vec,
                               alpha, beta, omega, n_alloc, distributions, N)
          ineq_jac <- if (nb > 1L) function(b_j) {
            Jc <- matrix(0, nrow = nb - 1L, ncol = nb)
            for (i in seq_len(nb - 1L)) { Jc[i, i] <- 1; Jc[i, i + 1L] <- -1 }
            Jc
          } else NULL
          tryCatch(
            nloptr::nloptr(x0 = b_j_init, eval_f = obj_fun,
                           eval_grad_f = grad_fun, lb = lb, ub = ub,
                           eval_g_ineq = ineq_fun, eval_jac_g_ineq = ineq_jac,
                           opts = list(algorithm = "NLOPT_LD_SLSQP",
                                       xtol_rel = 1e-8, ftol_rel = 1e-10,
                                       maxeval = 200L)),
            error = function(e) NULL)
        } else {
          tryCatch(
            nloptr::nloptr(x0 = b_j_init, eval_f = obj_fun, lb = lb, ub = ub,
                           eval_g_ineq = ineq_fun,
                           opts = list(algorithm = "NLOPT_LN_COBYLA",
                                       xtol_rel = 1e-6, ftol_rel = 1e-8,
                                       maxeval = 300L)),
            error = function(e) NULL)
        }

      if (!is.null(res_nlopt) && res_nlopt$objective < phi_before_j - 1e-14) {
        B[[j]] <- sort(res_nlopt$solution)
        sv     <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                              distributions, alpha, beta)
        W_h    <- sv$W_h
        N_h    <- sv$N_h
      }
    }

    # (B) Allocation update (Theorem 1(b)/(c), or CV-constrained programme)
    # Intermediate starts/boundary proposals may be CV-infeasible even when the
    # final selected design is feasible. Avoid surfacing warnings for those
    # transient candidates; feasibility remains attached to each allocation.
    n_alloc <- suppressWarnings(neyman_alloc(
      sigma2 = sv$sigma2, W_h = W_h, omega = omega, n_total = n_total,
      N_h = N_h, n_min = n_min, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit,
      cv_targets = cv_targets, mu = sv$mu))

    # Convergence (primary criterion: relative objective change)
    phi_new  <- .composite_obj(sv$sigma2, W_h, n_alloc, N_h, omega)
    phi_hist <- c(phi_hist, phi_new)
    if (abs(phi_new - phi_prev) / (1 + abs(phi_prev)) < tol_phi) {
      converged <- TRUE;  phi_prev <- phi_new;  break
    }
    phi_prev <- phi_new
  }

  list(B = B, n_alloc = n_alloc, phi = phi_prev, converged = converged,
       phi_trace = phi_hist)
}


# .phi_grad_boundary ----------------------------------------------------------
# Analytic gradient of Phi with respect to the interior boundaries of
# variable j.  The empirical Y-moments sigma2_hk AND the integer stratum
# sizes N_h are both piecewise constant in B (derivative zero almost
# everywhere), so the a.e.-exact gradient of the implemented objective
#   Phi = sum_h s_h W_h^2 (1/n_h - 1/N_h),  s_h = sum_k omega_k sigma2_hk
# is (matching eq. kkt_b of the paper with the sigma2 and N_h terms frozen):
#   dPhi/db = sum_h 2 s_h W_h (1/n_h - 1/N_h) dW_h/db
# with dW_h/db_l^{(j)} = +- f_j(b_l) * prod_{j' != j} DeltaF_{j'}  (Section 4.2 of the paper):
# the sign is + for the stratum with h_j = l (b_l is its upper endpoint)
# and - for the stratum with h_j = l + 1 (b_l is its lower endpoint).
.phi_grad_boundary <- function(b_j, j, B, data, study_vars, H_vec,
                               alpha, beta, omega, n_alloc, distributions, N) {
  b_j <- sort(b_j)
  B_cur <- B;  B_cur[[j]] <- b_j
  p  <- length(H_vec);  nb <- length(b_j)

  sv  <- .compute_stratum_means_vars(data, study_vars, B_cur, H_vec,
                                     distributions, alpha, beta)
  s_h <- as.numeric(sv$sigma2 %*% omega)       # sum_k omega_k sigma2_hk
  W_h <- sv$W_h

  # CDF increments per variable and per marginal stratum
  cdf_fns <- lapply(seq_len(p), function(jj)
    .build_cdf(distributions[[jj]], alpha[jj], beta[jj]))
  pdf_j   <- .build_pdf(distributions[[j]], alpha[j], beta[j])

  full_B <- vector("list", p)
  for (jj in seq_len(p)) full_B[[jj]] <- c(alpha[jj], B_cur[[jj]], beta[jj])
  incr <- vector("list", p)                    # DeltaF_{jj}(m), m = 1..H_jj
  for (jj in seq_len(p)) {
    Fv <- cdf_fns[[jj]](full_B[[jj]])
    incr[[jj]] <- pmax(diff(Fv), 0)
  }
  # Total mass (normalising constant), independent of interior boundaries
  s_norm <- prod(vapply(incr, sum, numeric(1)))
  if (s_norm < 1e-12) s_norm <- 1

  grids  <- lapply(seq_len(p), function(jj) seq_len(H_vec[jj]))
  combos <- as.matrix(expand.grid(grids, KEEP.OUT.ATTRS = FALSE))

  fpc_h  <- pmax(1 / pmax(n_alloc, 1) - 1 / pmax(sv$N_h, 1), 0)
  coef_h <- 2 * s_h * W_h * fpc_h

  grad <- numeric(nb)
  for (l in seq_len(nb)) {
    f_bl <- pdf_j(b_j[l])
    idx_up <- which(combos[, j] == l)          # b_l is upper endpoint: +
    idx_lo <- which(combos[, j] == l + 1L)     # b_l is lower endpoint: -
    g <- 0
    for (h in idx_up) {
      rest <- prod(vapply(seq_len(p)[-j],
                          function(jj) incr[[jj]][combos[h, jj]], numeric(1)))
      g <- g + coef_h[h] * ( f_bl * rest / s_norm)
    }
    for (h in idx_lo) {
      rest <- prod(vapply(seq_len(p)[-j],
                          function(jj) incr[[jj]][combos[h, jj]], numeric(1)))
      g <- g + coef_h[h] * (-f_bl * rest / s_norm)
    }
    grad[l] <- g
  }
  grad
}


# .kkt_residual ---------------------------------------------------------------
# Scaled KKT residual at the final solution (secondary diagnostic, paper
# Algorithm 1):
#   * allocation block: relative spread of the stationarity ratios
#     lambda_h = A_h^2 / n_h^2 over strata strictly inside their bounds
#     (all equal at an exact Neyman point);
#   * boundary block: max over interior boundaries of the scaled central
#     finite difference |dPhi/db| * (beta_j - alpha_j) / (1 + |Phi|).
.kkt_residual <- function(data, study_vars, B, n_alloc, H_vec, alpha, beta,
                          omega, distributions, n_min) {
  sv  <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                     distributions, alpha, beta)
  A_h <- .compute_A_h(sv$sigma2, sv$W_h, omega)
  n_a <- as.numeric(n_alloc)

  free <- n_a > n_min & n_a < sv$N_h
  r_alloc <- if (sum(free) >= 2L) {
    lam <- A_h[free]^2 / n_a[free]^2
    (max(lam) - min(lam)) / max(max(lam), 1e-300)
  } else 0

  phi0 <- .composite_obj(sv$sigma2, sv$W_h, n_a, sv$N_h, omega)
  r_bound <- 0
  for (j in seq_along(H_vec)) {
    db <- 1e-5 * (beta[j] - alpha[j])
    for (l in seq_along(B[[j]])) {
      for (sgn in c(-1, 1)) {
        Bt <- B;  Bt[[j]][l] <- Bt[[j]][l] + sgn * db
        if (sgn < 0) phi_m <- .composite_obj_at(Bt, data, study_vars, H_vec,
                                                distributions, alpha, beta,
                                                n_a, omega)
        else         phi_p <- .composite_obj_at(Bt, data, study_vars, H_vec,
                                                distributions, alpha, beta,
                                                n_a, omega)
      }
      d <- abs(phi_p - phi_m) / (2 * db) * (beta[j] - alpha[j]) / (1 + abs(phi0))
      r_bound <- max(r_bound, d)
    }
  }
  list(alloc = r_alloc, boundary = r_bound, max = max(r_alloc, r_bound))
}

.composite_obj_at <- function(B, data, study_vars, H_vec, distributions,
                              alpha, beta, n_alloc, omega) {
  sv <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                    distributions, alpha, beta)
  .composite_obj(sv$sigma2, sv$W_h, n_alloc, sv$N_h, omega)
}


.compute_stratum_stats <- function(data, study_vars, B, n_alloc,
                                   H_vec, omega, N, x_names, y_names,
                                   distributions = NULL,
                                   alpha = NULL, beta = NULL) {
  sv     <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                        distributions, alpha, beta)
  q      <- ncol(study_vars)
  labels <- .stratum_index_labels(H_vec)
  W_h    <- sv$W_h

  df <- data.frame(stratum = labels, N_h = sv$N_h, W_h = W_h,
                   n_h = as.integer(round(n_alloc)), stringsAsFactors = FALSE)
  for (k in seq_len(q)) {
    df[[paste0("mu_k",     k)]] <- sv$mu[, k]
    df[[paste0("sigma2_k", k)]] <- sv$sigma2[, k]
  }
  # FIX 1: corrected A_h (was the misaligned sweep with omega * W_h^2)
  df[["A_h"]] <- .compute_A_h(sv$sigma2, W_h, omega)
  attr(df, "strat_id") <- sv$strat_id
  df
}

# FIX 7: delegate to variance_lower_bound() so the reported gap and the
# exported bound share one implementation.
.optimality_gap <- function(strat_stats, phi_best, n_total, omega) {
  q          <- sum(grepl("^sigma2_k", names(strat_stats)))
  sigma2_mat <- as.matrix(
    strat_stats[, paste0("sigma2_k", seq_len(q)), drop = FALSE])
  lb <- variance_lower_bound(sigma2_mat, strat_stats$W_h, strat_stats$N_h,
                             omega, n_total)
  if (lb < 1e-15) return(0)
  (phi_best - lb) / abs(lb)
}
