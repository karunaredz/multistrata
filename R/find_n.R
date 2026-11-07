# =============================================================================
# find_n.R
#
# Inverse design: the smallest total sample size n that attains a set of
# coefficient-of-variation targets on every study variable, for a given
# stratification.  This is the operational question National Statistical
# Offices ask ("what is the cheapest survey that reaches CV <= 3% on all
# key indicators?") and is the dual of the fixed-n problem solved by
# stratify_multi().
#
# The routine reuses the CV-constrained allocation machinery already present
# in neyman_alloc() (constraint (C3) of Reddy & Khan, forthcoming): for fixed
# boundaries the per-stratum moments sigma2_hk, means mu_hk, weights W_h and
# sizes N_h do NOT depend on n, so they are computed once; feasibility of the
# CV targets is monotone in n, so the smallest feasible n is located by
# integer bisection.
#
# Reference:
#   Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#   Allocation for Multivariate Surveys (upcoming methods paper).
# =============================================================================


#' Smallest Sample Size Attaining Coefficient-of-Variation Targets
#'
#' Inverse (precision-driven) design.  Given a stratification and a set of
#' maximum coefficients of variation \eqn{\gamma_k} for the study-variable
#' mean estimators, \code{find_n()} returns the smallest total sample size
#' \eqn{n} for which an allocation exists that satisfies
#' \eqn{\mathrm{CV}(\hat{\mu}_k) \le \gamma_k} for every study variable
#' \eqn{k}, together with that allocation.
#'
#' For fixed strata boundaries the within-stratum moments, weights and
#' population sizes are independent of \eqn{n}; only the allocation changes.
#' Feasibility of the CV targets is monotone in \eqn{n} (more sample can only
#' lower each variance), so the minimal feasible \eqn{n} is found by integer
#' bisection over \code{n_range}.  At each trial size the convex CV-constrained
#' allocation sub-problem (constraint (C3) of Reddy & Khan, forthcoming) is solved via
#' \code{\link{neyman_alloc}}.
#'
#' Boundaries may be supplied through a pre-computed \code{"multistrata"}
#' object (\code{result}); otherwise \code{find_n()} first runs
#' \code{\link{stratify_multi}} once (at the upper end of \code{n_range}) to
#' obtain them.  With \code{reoptimize = TRUE} the boundaries are re-optimised
#' at every trial size (slower, occasionally slightly smaller \eqn{n}).
#'
#' @param data A \code{data.frame} or \code{matrix} of stratification
#'   variables (\eqn{N \times p}).
#' @param study_vars A \code{data.frame} or \code{matrix} of study variables
#'   (\eqn{N \times q}).
#' @param cv_targets Numeric vector of length \eqn{q}: the maximum acceptable
#'   coefficient of variation for each study-variable mean estimator (e.g.
#'   \code{c(0.03, 0.03)} for 3\%).
#' @param H_vec Integer vector of strata counts per stratification variable.
#'   Required unless \code{result} is supplied (then taken from it).
#' @param omega Study-variable weight vector (length \eqn{q}, summing to one).
#'   Defaults to equal weights (or the weights stored in \code{result}).
#' @param result Optional \code{"multistrata"} object supplying the strata
#'   boundaries (and \code{H_vec}, \code{omega}, \code{distributions}).  When
#'   \code{NULL}, boundaries are optimised internally.
#' @param distributions Parametric-mode specification passed through to
#'   \code{\link{stratify_multi}} / \code{\link{compute_stratum_stats_parametric}}
#'   (see \code{\link{stratify_multi}}).  Ignored when \code{result} is given
#'   (its own specification is reused).  Default \code{NULL} (empirical mode).
#' @param n_min Minimum sample size per stratum.  Default \code{2}.
#' @param n_range Optional length-2 integer vector \code{c(n_lo, n_hi)} giving
#'   the search bracket.  Defaults to \code{c(H_total * n_min, N)} where
#'   \eqn{N} is the population size (a full census is the precision ceiling).
#' @param reoptimize Logical.  Re-optimise boundaries at every trial \eqn{n}
#'   (via \code{stratify_multi(cv_targets = ...)}).  Default \code{FALSE}.
#' @param n_starts,max_iter Passed to \code{\link{stratify_multi}} when
#'   boundaries are optimised internally.  Defaults \code{20} and \code{200}.
#' @param verbose Print each bisection probe.  Default \code{FALSE}.
#'
#' @return An object of class \code{"findN"}: a list with components
#'   \code{n} (smallest feasible total sample size, or \code{NA} if the
#'   targets are unattainable within \code{n_range}), \code{allocation}
#'   (integer stratum sample sizes at \code{n}), \code{cv_achieved},
#'   \code{cv_targets}, \code{feasible}, \code{composite_variance},
#'   \code{boundaries}, \code{H_vec}, \code{omega}, \code{n_range}, and
#'   \code{trace} (a data.frame of the probed sizes and their feasibility).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "independent")
#' fn  <- find_n(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
#'               cv_targets = c(0.05, 0.05), H_vec = c(3, 3), n_starts = 5)
#' print(fn)
#' }
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{neyman_alloc}},
#'   \code{\link{design_cv}}
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#' Allocation for Multivariate Surveys (upcoming methods paper).
#'
#' @export
find_n <- function(data,
                   study_vars,
                   cv_targets,
                   H_vec         = NULL,
                   omega         = NULL,
                   result        = NULL,
                   distributions = NULL,
                   n_min         = 2L,
                   n_range       = NULL,
                   reoptimize    = FALSE,
                   n_starts      = 20L,
                   max_iter      = 200L,
                   verbose       = FALSE) {

  data       <- as.matrix(data)
  study_vars <- as.matrix(study_vars)
  N <- nrow(data);  q <- ncol(study_vars)

  if (missing(cv_targets) || is.null(cv_targets))
    stop("'cv_targets' must be supplied (one maximum CV per study variable).")
  if (length(cv_targets) != q)
    stop("'cv_targets' must have length equal to ncol(study_vars) = ", q, ".")
  if (any(cv_targets <= 0))
    stop("All elements of 'cv_targets' must be strictly positive.")

  # ------------------------------------------------------------------
  # 1.  Resolve boundaries / design metadata
  # ------------------------------------------------------------------
  if (!is.null(result)) {
    if (!inherits(result, "multistrata"))
      stop("'result' must be a 'multistrata' object from stratify_multi().")
    m             <- result$.meta
    H_vec         <- m$H_vec
    if (is.null(omega))         omega         <- m$omega
    if (is.null(distributions)) distributions <- m$distributions
    B             <- result$boundaries
    alpha         <- m$alpha
    beta          <- m$beta
  } else {
    if (is.null(H_vec))
      stop("'H_vec' must be supplied when 'result' is NULL.")
    if (length(H_vec) != ncol(data))
      stop("'H_vec' must have one entry per column of 'data'.")
    if (is.null(omega)) omega <- rep(1 / q, q)
    alpha <- apply(data, 2L, min)
    beta  <- apply(data, 2L, max)
  }
  if (abs(sum(omega) - 1) > 1e-8) stop("'omega' must sum to 1.")

  H_total <- prod(H_vec)

  # Default search bracket: [H_total * n_min, N]  (census = precision ceiling)
  if (is.null(n_range)) n_range <- c(H_total * n_min, N)
  n_range <- as.integer(round(n_range))
  n_lo <- max(n_range[1L], H_total * n_min)
  n_hi <- min(n_range[2L], N)
  if (n_hi <= n_lo)
    stop("Empty search bracket: n_range = (", n_lo, ", ", n_hi, ").")

  B <- if (!is.null(result)) result$boundaries else NULL

  # If we must optimise boundaries and none supplied, do it once at n_hi.
  if (is.null(result) && !reoptimize) {
    if (verbose) message("  [find_n] optimising boundaries once at n = ", n_hi, " ...")
    result <- stratify_multi(
      data = data, study_vars = study_vars, H_vec = H_vec, n_total = n_hi,
      omega = omega, n_min = n_min, n_starts = n_starts, max_iter = max_iter,
      distributions = distributions, verbose = FALSE)
    B <- result$boundaries
  }

  # ------------------------------------------------------------------
  # 2.  Feasibility oracle at a trial size n
  # ------------------------------------------------------------------
  if (!reoptimize) {
    # Fixed boundaries: moments computed once, invariant to n
    sv     <- .compute_stratum_means_vars(data, study_vars, B, H_vec,
                                          distributions, alpha, beta)
    sigma2 <- sv$sigma2;  W_h <- sv$W_h;  N_h <- sv$N_h;  mu <- sv$mu

    feasible_at <- function(n) {
      alloc <- suppressWarnings(neyman_alloc(
        sigma2 = sigma2, W_h = W_h, omega = omega, n_total = n, N_h = N_h,
        n_min = n_min, cv_targets = cv_targets, mu = mu))
      cv  <- attr(alloc, "cv_achieved")
      ok  <- isTRUE(attr(alloc, "cv_feasible")) &&
             all(cv <= cv_targets + 1e-6)
      list(feasible = ok, alloc = alloc, cv = cv)
    }
  } else {
    feasible_at <- function(n) {
      res <- suppressWarnings(stratify_multi(
        data = data, study_vars = study_vars, H_vec = H_vec, n_total = n,
        omega = omega, n_min = n_min, n_starts = n_starts, max_iter = max_iter,
        cv_targets = cv_targets, distributions = distributions, verbose = FALSE))
      cv <- res$cv_achieved
      ok <- !is.null(cv) && all(cv <= cv_targets + 1e-6)
      list(feasible = ok, alloc = res$allocation, cv = cv, result = res)
    }
  }

  trace <- data.frame(n = integer(0), feasible = logical(0))
  probe <- function(n) {
    r <- feasible_at(n)
    trace[nrow(trace) + 1L, ] <<- list(as.integer(n), r$feasible)
    if (verbose)
      message(sprintf("  [find_n] n = %6d  ->  %s  (CV = %s)",
                      as.integer(n),
                      if (r$feasible) "feasible" else "infeasible",
                      if (is.null(r$cv)) "NA"
                      else paste(sprintf("%.4f", r$cv), collapse = ", ")))
    r
  }

  # ------------------------------------------------------------------
  # 3.  Bracket + integer bisection for the smallest feasible n
  # ------------------------------------------------------------------
  r_hi <- probe(n_hi)
  if (!r_hi$feasible) {
    warning("CV targets are unattainable within n_range (infeasible even at ",
            "n = ", n_hi, ", the upper bracket). Returning n = NA.")
    return(.make_findN(NA_integer_, r_hi$alloc, r_hi$cv, cv_targets, FALSE,
                       omega = omega, B = B, H_vec = H_vec,
                       n_range = c(n_lo, n_hi),
                       trace = trace[order(trace$n), , drop = FALSE]))
  }

  r_lo <- probe(n_lo)
  if (r_lo$feasible) {
    best_n <- n_lo;  best_r <- r_lo
  } else {
    lo <- n_lo;  hi <- n_hi
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      rm  <- probe(mid)
      if (rm$feasible) hi <- mid else lo <- mid
    }
    best_n <- hi;  best_r <- probe(hi)
  }

  # Composite variance at the solution
  if (!reoptimize) {
    phi <- composite_variance(sigma2, W_h, best_r$alloc, N_h, omega)
    B_out <- B
  } else {
    phi   <- best_r$result$objective
    B_out <- best_r$result$boundaries
  }

  .make_findN(as.integer(best_n), best_r$alloc, best_r$cv, cv_targets, TRUE,
              omega = omega, B = B_out, H_vec = H_vec, phi = phi,
              n_range = c(n_lo, n_hi),
              trace = trace[order(trace$n), , drop = FALSE])
}


# .make_findN -----------------------------------------------------------------
.make_findN <- function(n, alloc, cv, cv_targets, feasible,
                        omega = NULL, B = NULL, H_vec = NULL, phi = NA_real_,
                        n_range = NULL, trace = NULL) {
  structure(
    list(
      n                  = n,
      allocation         = alloc,
      cv_achieved        = cv,
      cv_targets         = cv_targets,
      feasible           = feasible,
      composite_variance = phi,
      boundaries         = B,
      H_vec              = H_vec,
      omega              = omega,
      n_range            = n_range,
      trace              = trace
    ),
    class = "findN"
  )
}


#' Print method for findN
#' @param x A \code{"findN"} object.
#' @param ... Unused.
#' @export
print.findN <- function(x, ...) {
  cat("\n")
  cat("== find_n(): precision-driven sample size ======\n")
  cat(sprintf("  Strata per variable : %s  (H_total = %d)\n",
              paste(x$H_vec, collapse = " x "), prod(x$H_vec)))
  cat(sprintf("  CV targets          : %s\n",
              paste(sprintf("%.2f%%", 100 * x$cv_targets), collapse = ", ")))
  cat(sprintf("  Search bracket      : [%d, %d]\n",
              x$n_range[1L], x$n_range[2L]))
  cat("------------------------------------------------\n")
  if (isTRUE(x$feasible)) {
    cat(sprintf("  Smallest feasible n : %d\n", x$n))
    cat(sprintf("  CV achieved         : %s\n",
                paste(sprintf("%.2f%%", 100 * x$cv_achieved), collapse = ", ")))
    if (is.finite(x$composite_variance))
      cat(sprintf("  Composite variance  : %.8f\n", x$composite_variance))
    cat(sprintf("  Probes evaluated    : %d\n", nrow(x$trace)))
  } else {
    cat("  CV targets UNATTAINABLE within the search bracket.\n")
    cat("  Consider a finer stratification (larger H_vec) or relaxed targets.\n")
  }
  cat("================================================\n\n")
  invisible(x)
}
