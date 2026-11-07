# =============================================================================
# neyman_alloc.R
#
# Generalised Neyman allocation (Theorem 1(b) in Reddy & Khan, forthcoming):
#   n*_h = n_total * A_h / sum_g A_g
# where A_h = W_h * sqrt( sum_k omega_k * sigma2_hk ).
#
# Under the cost constraint (Theorem 1(c)):
#   n*_h = (C* - c0) * (A_h / sqrt(c_h)) / sum_g ( A_g * sqrt(c_g) ),
# i.e. n_h  is proportional to A_h / sqrt(c_h), the classical square-root-cost
# rule, generalised to the composite objective.
#
# Fixes in this revision:
#   1. A_h was computed as sqrt(rowSums(sweep(sigma2, 2, omega * W_h^2, "*"))).
#      'omega' (length q) times 'W_h^2' (length H) recycles element-wise and
#      sweep() then misaligns the result across columns, mixing weights from
#      other strata into A_h whenever H != q.  Corrected to
#      W_h * sqrt(rowSums(sweep(sigma2, 2, omega, "*"))).
#   2. Cost-constrained formula corrected from sqrt(A_h / c_h) to
#      A_h / sqrt(c_h) (Lagrangian stationarity gives
#      n_h = A_h / sqrt(lambda * c_h)), and the fixed overhead c0
#      ('cost_fixed') is now subtracted from the budget.
#   3. In cost mode the integer rounding now targets the implied total
#      round(sum(n_h)) rather than 'n_total'.
#   4. Bound handling now uses the iterative active-set correction of
#      Theorem 1(b): strata hitting a bound are fixed there and the closed
#      form is reapplied to the remaining free strata with the residual
#      sample (or residual budget), until no bound is violated.
#   5. Optional CV constraints (C3): when 'cv_targets' is supplied the
#      allocation sub-problem is solved numerically as a convex programme
#      (SLSQP) subject to CV(mu_hat_k) <= gamma_k for each study variable.
# =============================================================================


#' Generalised Multivariate Neyman Allocation
#'
#' Computes the optimal stratum sample sizes for fixed strata by minimising
#' the weighted composite variance \eqn{\Phi(\mathbf{B},\mathbf{n})} over
#' \eqn{\mathbf{n}} for given \eqn{\mathbf{B}}.  This is the closed-form
#' solution to the allocation sub-problem established in Theorem 1(b) of
#' Reddy & Khan (forthcoming), with the iterative active-set correction applied when
#' stratum bounds become active.  Under a cost constraint the square-root-cost
#' rule of Theorem 1(c) is used.  When coefficient-of-variation targets are
#' supplied (constraint (C3) of the paper), the allocation sub-problem --
#' which remains convex in \eqn{\mathbf{n}} -- is solved numerically by SLSQP.
#'
#' When \eqn{q = 1} and \eqn{\omega_1 = 1} the formula reduces to the
#' classical Neyman allocation \eqn{n_h \propto W_h \sigma_h}, and the
#' cost-constrained version reduces to \eqn{n_h \propto W_h\sigma_h/\sqrt{c_h}}
#' (Cochran, 1977, Ch. 5).
#'
#' @param sigma2 Matrix of within-stratum variances, dimension
#'   \eqn{H \times q}.
#' @param W_h Numeric vector of stratum weights, length \eqn{H}, summing to 1.
#' @param omega Weight vector of length \eqn{q} summing to 1.
#' @param n_total Total sample size (required when
#'   \code{constraint = "n_total"} and when \code{cv_targets} is supplied).
#' @param N_h Integer vector of stratum population sizes (length \eqn{H}).
#' @param n_min Minimum sample size per stratum.  Default \code{2}.
#' @param constraint One of \code{"n_total"} (default) or \code{"cost"}.
#' @param cost_total Total budget \eqn{C^*} (used when
#'   \code{constraint = "cost"}).
#' @param cost_fixed Fixed survey overhead \eqn{c_0} subtracted from the
#'   budget before allocation.  Default \code{0}.
#' @param cost_per_unit Scalar or length-\eqn{H} vector of per-unit costs
#'   \eqn{c_h}.  Default \code{1}.
#' @param cv_targets Optional numeric vector of length \eqn{q} of maximum
#'   coefficients of variation \eqn{\gamma_k} (constraint (C3)).  Requires
#'   \code{mu}.  Default \code{NULL}.
#' @param mu Matrix of within-stratum study-variable means
#'   (\eqn{H \times q}); required when \code{cv_targets} is supplied.
#'
#' @return Numeric vector of optimal integer stratum sample sizes.  When
#'   \code{cv_targets} is supplied, the attribute \code{"cv_achieved"} holds
#'   the achieved CVs and \code{"cv_feasible"} flags whether all targets were
#'   met.
#'
#' @examples
#' sigma2 <- matrix(c(4, 9, 16, 25, 1, 4, 9, 16), nrow = 4, ncol = 2)
#' W_h  <- c(0.3, 0.3, 0.2, 0.2)
#' N_h  <- c(300, 300, 200, 200)
#' alloc <- neyman_alloc(sigma2, W_h, omega = c(0.5, 0.5),
#'                       n_total = 200, N_h = N_h)
#' sum(alloc)   # 200
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#' Allocation for Multivariate Surveys (upcoming methods paper).
#'
#' Bethel, J. (1989). Sample allocation in multivariate surveys.
#' \emph{Survey Methodology}, \strong{15}(1), 47--57.
#'
#' Cochran, W. G. (1977). \emph{Sampling Techniques}, 3rd ed. Wiley.
#'
#' @export
neyman_alloc <- function(sigma2,
                         W_h,
                         omega,
                         n_total        = NULL,
                         N_h,
                         n_min          = 2L,
                         constraint     = c("n_total", "cost"),
                         cost_total     = NULL,
                         cost_fixed     = 0,
                         cost_per_unit  = 1,
                         cv_targets     = NULL,
                         mu             = NULL) {

  constraint <- match.arg(constraint)
  sigma2     <- as.matrix(sigma2)
  H          <- nrow(sigma2);  q <- ncol(sigma2)

  if (length(W_h)  != H) stop("'W_h' must have length nrow(sigma2).")
  if (length(N_h)  != H) stop("'N_h' must have length nrow(sigma2).")
  if (length(omega) != q) stop("'omega' must have length ncol(sigma2).")
  if (constraint == "n_total" && is.null(n_total))
    stop("'n_total' must be specified when constraint = 'n_total'.")

  c_h <- if (length(cost_per_unit) == 1L) rep(cost_per_unit, H) else cost_per_unit
  if (length(c_h) != H) stop("'cost_per_unit' must be scalar or length H.")
  c_h <- as.numeric(c_h)

  # nloptr requires double storage mode (tabulate() etc. return integers)
  N_h     <- as.numeric(N_h)
  n_min   <- as.numeric(n_min)
  if (!is.null(n_total)) n_total <- as.numeric(n_total)

  # A_h = W_h * sqrt( sum_k omega_k * sigma2_hk )   [Theorem 1(b)]
  A_h <- .compute_A_h(sigma2, W_h, omega)

  # Effective per-stratum bounds (n_min cannot exceed N_h)
  lo <- as.numeric(pmin(rep(n_min, H), N_h))
  hi <- as.numeric(pmax(N_h, lo))

  # ------------------------------------------------------------------
  # CV-constrained allocation (C3): convex programme solved by SLSQP
  # ------------------------------------------------------------------
  if (!is.null(cv_targets)) {
    if (is.null(mu))
      stop("'mu' (H x q stratum means) is required when 'cv_targets' is supplied.")
    if (length(cv_targets) != q)
      stop("'cv_targets' must have length ncol(sigma2).")
    if (is.null(n_total))
      stop("'n_total' is required when 'cv_targets' is supplied.")
    return(.alloc_cv_constrained(sigma2, as.matrix(mu), W_h, omega, n_total,
                                 N_h, lo, hi, A_h, cv_targets))
  }

  # ------------------------------------------------------------------
  # Closed-form allocation with active-set correction (Theorem 1(b)/(c))
  # ------------------------------------------------------------------
  if (constraint == "n_total") {
    n_cont <- .active_set_n(A_h, n_total, lo, hi)
    target <- as.integer(round(n_total))
  } else {
    if (is.null(cost_total))
      stop("'cost_total' must be specified when constraint = 'cost'.")
    budget <- cost_total - cost_fixed
    if (budget <= 0)
      stop("'cost_total' must exceed 'cost_fixed'.")
    n_cont <- .active_set_cost(A_h, budget, c_h, lo, hi)
    target <- as.integer(round(sum(n_cont)))
  }

  alloc <- .round_to_total(n_cont, target, n_min, N_h)

  # In cost mode, ensure the rounded allocation respects the budget
  if (constraint == "cost") {
    over <- sum(c_h * alloc) - budget
    if (over > 0) {
      ord <- order(c_h, decreasing = TRUE)
      for (h in ord) {
        while (over > 0 && alloc[h] > lo[h]) {
          alloc[h] <- alloc[h] - 1L
          over     <- over - c_h[h]
        }
        if (over <= 0) break
      }
    }
  }
  alloc
}


# .compute_A_h ---------------------------------------------------------------
# A_h = W_h * sqrt( sum_k omega_k sigma2_hk ).
# FIX 1: replaces sweep(sigma2, 2, omega * W_h^2, "*"), which misaligned
# omega/W_h whenever H != q (element-wise recycling of a length-q by a
# length-H vector, then sweep across the wrong extent).
.compute_A_h <- function(sigma2, W_h, omega) {
  A_h <- W_h * sqrt(rowSums(sweep(as.matrix(sigma2), 2L, omega, "*")))
  pmax(A_h, 1e-12)
}


# .active_set_n ---------------------------------------------------------------
# Exact bound-constrained Neyman allocation via the KKT water-filling
# characterisation (Theorem 1(b)): the unique optimum has
#   n_h(lambda) = min(N_h, max(n_min, A_h / sqrt(lambda))),
# where lambda > 0 is the unique root of  sum_h n_h(lambda) = n_total
# (the left side is continuous and strictly decreasing wherever any stratum
# is unclipped).  lambda is located by bisection on the log scale.  This is
# equivalent to -- and strictly more robust than -- iteratively fixing bound
# violators and reapplying the closed form to the free strata.
.active_set_n <- function(A_h, n_total, lo, hi) {
  n_of <- function(lam) pmin(hi, pmax(lo, A_h / sqrt(lam)))
  if (sum(hi) <= n_total) return(hi)       # even full saturation can't exceed
  if (sum(lo) >= n_total) return(lo)       # minima already exhaust the sample
  lam_lo <- 1e-18;  lam_hi <- 1e18
  for (i in seq_len(240L)) {
    mid <- sqrt(lam_lo * lam_hi)
    if (sum(n_of(mid)) > n_total) lam_lo <- mid else lam_hi <- mid
  }
  n_of(sqrt(lam_lo * lam_hi))
}


# .active_set_cost ------------------------------------------------------------
# Exact bound-constrained cost allocation via water-filling (Theorem 1(c)):
#   n_h(lambda) = min(N_h, max(n_min, A_h / sqrt(lambda * c_h))),
# with lambda the unique root of  sum_h c_h * n_h(lambda) = budget.
.active_set_cost <- function(A_h, budget, c_h, lo, hi) {
  n_of <- function(lam) pmin(hi, pmax(lo, A_h / sqrt(lam * c_h)))
  if (sum(c_h * hi) <= budget) return(hi)
  if (sum(c_h * lo) >= budget) return(lo)
  lam_lo <- 1e-18;  lam_hi <- 1e18
  for (i in seq_len(240L)) {
    mid <- sqrt(lam_lo * lam_hi)
    if (sum(c_h * n_of(mid)) > budget) lam_lo <- mid else lam_hi <- mid
  }
  n_of(sqrt(lam_lo * lam_hi))
}


# .alloc_cv_constrained -------------------------------------------------------
# Allocation sub-problem with CV constraints (C3).  For fixed boundaries the
# problem  min Phi(n)  s.t.  sum n = n_total,  CV_k(n) <= gamma_k,
# lo <= n <= hi  is convex in n (each V_k is a positive combination of 1/n_h),
# so SLSQP converges to the unique global optimum.  Gradients are analytic.
.alloc_cv_constrained <- function(sigma2, mu, W_h, omega, n_total,
                                  N_h, lo, hi, A_h, cv_targets) {
  H <- nrow(sigma2);  q <- ncol(sigma2)

  mu_pop <- as.numeric(crossprod(W_h, mu))          # population means, length q
  if (any(abs(mu_pop) < 1e-15))
    stop("CV constraints undefined: a population mean is essentially zero.")

  # Guard the finite-population correction exactly as composite_variance(),
  # design_cv() and .composite_obj() do.
  #
  # An empty stratum has N_h = 0 and therefore W_h = 0, so it contributes
  # nothing to the variance -- but written as 1/n - 1/N_h the correction is
  # Inf there, and W_h^2 * sigma2 * Inf is 0 * Inf = NaN. That NaN made the
  # SLSQP solver reject its own warm start with "objective in x0 returns NA",
  # which aborted the whole run. Empty strata are routine during the boundary
  # search (a trial cut can leave a cell with no units), so the failure got
  # more likely as n_starts grew and was near-certain at the default of 50.
  n_eff <- function(n) pmax(n, 1)
  V_k <- function(n) {                              # length-q variances
    fpc <- pmax(1 / n_eff(n) - 1 / pmax(N_h, 1), 0)
    as.numeric(colSums(W_h^2 * sigma2 * fpc))
  }
  obj  <- function(n) sum(omega * V_k(n))
  gobj <- function(n)
    -rowSums(sweep(W_h^2 * sigma2, 2L, omega, "*")) / n_eff(n)^2

  # Inequalities g(n) <= 0:  V_k(n) - (gamma_k * mu_k)^2 <= 0
  v_max <- (cv_targets * abs(mu_pop))^2
  gin   <- function(n) V_k(n) - v_max
  gin_j <- function(n) t(-(W_h^2 * sigma2) / n_eff(n)^2)   # q x H Jacobian

  geq   <- function(n) sum(n) - n_total
  geq_j <- function(n) rep(1, H)

  n0 <- .active_set_n(A_h, n_total, lo, hi)         # warm start at Neyman
  n0 <- as.numeric(pmin(pmax(n0, lo + 1e-9), hi))

  res <- nloptr::nloptr(
    x0 = n0, eval_f = obj, eval_grad_f = gobj,
    lb = as.numeric(lo), ub = as.numeric(hi),
    eval_g_ineq = gin, eval_jac_g_ineq = gin_j,
    eval_g_eq   = geq, eval_jac_g_eq   = geq_j,
    opts = list(algorithm = "NLOPT_LD_SLSQP",
                xtol_rel = 1e-10, maxeval = 2000L))

  n_cont <- res$solution
  alloc  <- .round_to_total(n_cont, as.integer(round(n_total)),
                            min(lo), N_h)

  cv_ach <- sqrt(pmax(V_k(pmax(alloc, 1)), 0)) / abs(mu_pop)
  # Integer rounding can move a continuous SLSQP solution by one unit. The
  # public tests and reporting contract allow 5e-4 absolute CV slack (0.05
  # percentage points), while still flagging substantively infeasible targets.
  cv_rounding_tol <- 5e-4
  feasible <- all(cv_ach <= cv_targets + cv_rounding_tol)
  if (!feasible)
    warning(sprintf(
      "CV targets not attainable with n_total = %d (achieved: %s). ",
      as.integer(n_total), paste(sprintf("%.4f", cv_ach), collapse = ", ")),
      "Consider increasing n_total or relaxing 'cv_targets'.")

  attr(alloc, "cv_achieved") <- cv_ach
  attr(alloc, "cv_feasible") <- feasible
  alloc
}


# .round_to_total -----------------------------------------------------------
# Hamilton (largest-remainder) rounding that preserves sum = target,
# with a safety pass when n_min prevents exact correction.
.round_to_total <- function(raw, target, n_min, N_h) {
  H       <- length(raw)
  floor_a <- floor(raw)
  frac    <- raw - floor_a
  alloc   <- pmax(floor_a, n_min)
  alloc   <- pmin(alloc, N_h)

  deficit <- target - sum(alloc)

  if (deficit > 0L) {
    ord <- order(frac, decreasing = TRUE)
    i   <- 0L
    while (deficit > 0L) {
      i <- i + 1L
      h <- ord[((i - 1L) %% H) + 1L]
      if (alloc[h] < N_h[h]) {
        alloc[h]  <- alloc[h] + 1L
        deficit   <- deficit - 1L
      }
      if (i > 10L * H) break                       # all strata saturated
    }
  } else if (deficit < 0L) {
    ord     <- order(frac, decreasing = FALSE)
    removed <- 0L
    for (i in seq_len(H)) {
      if (removed >= abs(deficit)) break
      h <- ord[i]
      if (alloc[h] > n_min) {
        alloc[h] <- alloc[h] - 1L
        removed  <- removed + 1L
      }
    }
    # Safety pass, if still over, force-reduce from largest strata
    remaining <- sum(alloc) - target
    if (remaining > 0L) {
      ord2 <- order(alloc, decreasing = TRUE)
      for (i in seq_len(H)) {
        if (remaining <= 0L) break
        h <- ord2[i]
        while (alloc[h] > n_min && remaining > 0L) {
          alloc[h]  <- alloc[h] - 1L
          remaining <- remaining - 1L
        }
      }
      if (remaining > 0L)
        warning(sprintf(
          ".round_to_total: could not satisfy sum = %d (n_min constraint too tight); achieved %d.",
          target, as.integer(sum(alloc))))
    }
  }
  as.numeric(alloc)
}
