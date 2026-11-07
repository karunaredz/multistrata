# =============================================================================
# omega_sweep.R
#
# Data-mode counterpart of omega_sweep_density() (density_analysis.R).
#
# Every other analysis tool in the package ships as a data/density pair --
# find_n()/find_n_density(), choose_H()/choose_H_density(),
# benchmark_designs()/benchmark_density_designs() -- and the omega sweep was
# the one asymmetry: only the density branch was exported, so the empirical
# and parametric modes had no supported way to trace the priority-weight
# trade-off. This file closes that gap and is the single implementation used
# by both the console API and the bundled Shiny app's omega explorer.
# =============================================================================


#' Trade-Off Between Two Study Variables Across Priority Weights
#'
#' Re-optimises a design over a grid of priority weights for two study
#' variables and returns objective and precision diagnostics at each grid
#' point, so the analyst can see what protecting one variable's precision
#' costs the other.
#'
#' The weight assigned to \code{pair[1]} is swept across \code{grid}, expressed
#' as its share of the two variables' combined weight. The remaining
#' variables keep their weights, and the full vector is renormalised to sum to
#' one at every grid point, so only the \emph{relative} priority within the
#' pair changes. The design is re-optimised from scratch at each point:
#' boundaries move with the weights, which is precisely why the trade-off
#' cannot be read off a single fitted design.
#'
#' This is the empirical/parametric analogue of
#' \code{\link{omega_sweep_density}} and returns the same columns, so the two
#' are interchangeable downstream.
#'
#' @param data       Matrix/data.frame of stratification variables
#'   (\eqn{N \times p}).
#' @param study_vars Matrix/data.frame of study variables (\eqn{N \times q}).
#' @param H_vec      Integer vector of strata counts per stratification
#'   variable. May be \code{NULL} when \code{result} is supplied.
#' @param n_total    Total sample size.
#' @param omega      Baseline study-variable weight vector (length \eqn{q}).
#'   Defaults to equal weights, or to the weights of \code{result}.
#' @param pair       Integer vector of length 2, or a character vector of two
#'   study-variable names, giving the variables whose relative weights are
#'   swept. Defaults to the first two.
#' @param grid       Numeric vector of shares in \eqn{(0, 1)} assigned to
#'   \code{pair[1]} within the pair's combined weight.
#' @param constraint Either \code{"n_total"} or \code{"cost"}.
#' @param cost_total Total budget when \code{constraint = "cost"}.
#' @param cost_fixed Fixed survey overhead.
#' @param cost_per_unit Per-unit cost, scalar or one value per stratum.
#' @param cv_targets Optional coefficient-of-variation targets, one per study
#'   variable.
#' @param n_min      Minimum sample size per stratum.
#' @param n_starts   Multi-starts per grid point. Kept modest by default
#'   because the whole grid is re-optimised.
#' @param max_iter   Maximum alternating-optimisation iterations per start.
#' @param method     Boundary solver, \code{"cobyla"} or \code{"slsqp"}.
#' @param distributions Optional parametric marginal specifications, as
#'   accepted by \code{\link{stratify_multi}}.
#' @param result     Optional fitted \code{"multistrata"} object whose
#'   \code{H_vec}, \code{omega}, \code{method} and \code{distributions} are
#'   reused, so a sweep can follow a fitted design without restating it.
#' @param verbose    Logical; report progress per grid point.
#'
#' @return A \code{data.frame} with one row per grid point and columns
#' \describe{
#'   \item{\code{w}}{weight assigned to \code{pair[1]} after renormalisation.}
#'   \item{\code{phi}}{optimal composite variance \eqn{\Phi^*} at that weight.}
#'   \item{\code{cv1}, \code{cv2}}{achieved coefficients of variation (per
#'     cent) of \code{pair[1]} and \code{pair[2]}.}
#'   \item{\code{gap}}{relative optimality gap (per cent).}
#' }
#' The names of the swept variables are attached as the \code{"pair"} and
#' \code{"pair_names"} attributes.
#'
#' @seealso \code{\link{omega_sweep_density}} for the data-free analogue,
#'   \code{\link{stratify_multi}}, \code{\link{design_cv}}
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#' Allocation for Multivariate Surveys (upcoming methods paper).
#'
#' @examples
#' set.seed(99)
#' pop <- sim_population(N = 3000, p = 2, q = 2, scenario = "correlated")
#' sw <- omega_sweep(
#'   data       = pop[, c("X1", "X2")],
#'   study_vars = pop[, c("Y1", "Y2")],
#'   H_vec      = c(2, 2),
#'   n_total    = 200,
#'   grid       = c(0.25, 0.5, 0.75),
#'   n_starts   = 1
#' )
#' sw
#'
#' @export
omega_sweep <- function(data,
                        study_vars,
                        H_vec         = NULL,
                        n_total,
                        omega         = NULL,
                        pair          = c(1L, 2L),
                        grid          = seq(0.05, 0.95, length.out = 7L),
                        constraint    = c("n_total", "cost"),
                        cost_total    = NULL,
                        cost_fixed    = 0,
                        cost_per_unit = 1,
                        cv_targets    = NULL,
                        n_min         = 2L,
                        n_starts      = 5L,
                        max_iter      = 60L,
                        method        = c("cobyla", "slsqp"),
                        distributions = NULL,
                        result        = NULL,
                        verbose       = FALSE) {

  constraint <- match.arg(constraint)

  data       <- as.matrix(data)
  study_vars <- as.matrix(study_vars)
  if (nrow(study_vars) != nrow(data))
    stop("'data' and 'study_vars' must have the same number of rows.")
  q <- ncol(study_vars)
  y_names <- colnames(study_vars)
  if (is.null(y_names)) y_names <- paste0("Y", seq_len(q))

  # Inherit design metadata from a fitted object, as find_n() does.
  if (!is.null(result)) {
    if (!inherits(result, "multistrata"))
      stop("'result' must be a 'multistrata' object from stratify_multi().")
    m <- result$.meta
    if (is.null(H_vec))         H_vec         <- m$H_vec
    if (is.null(omega))         omega         <- m$omega
    if (is.null(distributions)) distributions <- m$distributions
    if (missing(method) && !is.null(m$method)) method <- m$method
  }
  method <- match.arg(method)

  if (is.null(H_vec))
    stop("'H_vec' must be supplied when 'result' is NULL.")
  if (length(H_vec) != ncol(data))
    stop("'H_vec' must have one entry per column of 'data'.")
  if (missing(n_total) || is.null(n_total))
    stop("'n_total' must be supplied.")

  if (is.null(omega)) omega <- rep(1 / q, q)
  if (length(omega) != q)
    stop("'omega' must have length equal to ncol(study_vars) = ", q, ".")

  # Accept either indices or study-variable names for the swept pair.
  if (is.character(pair)) {
    idx <- match(pair, y_names)
    if (anyNA(idx))
      stop("'pair' names not found among the study variables: ",
           paste(pair[is.na(idx)], collapse = ", "), ".")
    pair <- idx
  }
  pair <- as.integer(pair)
  if (length(pair) != 2L || anyNA(pair) ||
      any(pair < 1L) || any(pair > q) || pair[1L] == pair[2L])
    stop("'pair' must name or index two different study variables.")

  grid <- as.numeric(grid)
  if (!length(grid) || any(!is.finite(grid)) || any(grid <= 0) || any(grid >= 1))
    stop("'grid' must contain finite shares strictly between 0 and 1.")

  total_pair <- sum(omega[pair])
  if (!is.finite(total_pair) || total_pair <= 0)
    stop("The swept pair must carry strictly positive combined weight.")

  out <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    om <- omega
    om[pair[1L]] <- grid[i] * total_pair
    om[pair[2L]] <- (1 - grid[i]) * total_pair
    om <- om / sum(om)

    if (isTRUE(verbose))
      message(sprintf("  [omega_sweep] %d/%d: omega(%s) = %.3f",
                      i, length(grid), y_names[pair[1L]], om[pair[1L]]))

    r <- stratify_multi(
      data = data, study_vars = study_vars, H_vec = H_vec,
      n_total = n_total, omega = om, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit, cv_targets = cv_targets, n_min = n_min,
      n_starts = n_starts, max_iter = max_iter, method = method,
      distributions = distributions, verbose = FALSE)

    mats <- .extract_stratum_moments(r$stratum_stats, q)
    cv <- as.numeric(design_cv(mats$sigma2, mats$mu, r$stratum_stats$W_h,
                               r$allocation, r$stratum_stats$N_h))
    out[[i]] <- data.frame(
      w   = om[pair[1L]],
      phi = r$objective,
      cv1 = 100 * cv[pair[1L]],
      cv2 = 100 * cv[pair[2L]],
      gap = 100 * r$optimality_gap,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, out)
  attr(tab, "pair")       <- pair
  attr(tab, "pair_names") <- y_names[pair]
  tab
}


# Pull the per-stratum study-variable moment matrices out of a stratum_stats
# data frame. Shared by omega_sweep() and the density analysis tools; the
# column naming (sigma2_k1.., mu_k1..) is fixed by stratify_multi().
.extract_stratum_moments <- function(ss, q) {
  list(
    sigma2 = as.matrix(ss[, paste0("sigma2_k", seq_len(q)), drop = FALSE]),
    mu     = as.matrix(ss[, paste0("mu_k",     seq_len(q)), drop = FALSE])
  )
}