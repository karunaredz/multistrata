# =============================================================================
# choose_H.R
#
# Number-of-strata selector.  "How many strata?" is the practitioner's second
# question (after "which variables?").  choose_H() traces the optimal
# composite variance Phi* against a set of candidate strata grids and locates
# the elbow -- the point beyond which further subdivision buys negligible
# precision, so the extra field/operational cost of more strata is not
# justified.  This mirrors the diagnostic pattern proven useful in the
# stratifyR 2.0 Shiny application.
#
# Reference:
#   Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#   Allocation for Multivariate Surveys (upcoming methods paper).
# =============================================================================


#' Choose the Number of Strata by Tracing Phi*
#'
#' Fits \code{\link{stratify_multi}} across a family of candidate strata grids
#' and returns, for each, the optimal composite variance \eqn{\Phi^*}, the
#' optimality gap and the relative precision gain over the coarsest grid.
#' A recommended grid is flagged at the elbow of the \eqn{\Phi^*} curve -- the
#' smallest grid beyond which the marginal reduction in \eqn{\Phi^*} falls
#' below \code{tol_gain}.
#'
#' @param data A \code{data.frame} or \code{matrix} of stratification
#'   variables (\eqn{N \times p}).
#' @param study_vars A \code{data.frame} or \code{matrix} of study variables
#'   (\eqn{N \times q}).
#' @param n_total Total sample size.
#' @param H_list A list of candidate \code{H_vec} integer vectors (each of
#'   length \eqn{p}).  If \code{NULL}, a default symmetric ladder is built from
#'   \code{H_each}: \code{lapply(H_each, function(h) rep(h, p))}.
#' @param H_each Integer vector of per-variable strata counts used to build a
#'   symmetric ladder when \code{H_list} is \code{NULL}.  Default
#'   \code{2:5}.
#' @param omega Study-variable weight vector (length \eqn{q}).  Defaults to
#'   equal weights.
#' @param distributions Parametric-mode specification passed to
#'   \code{\link{stratify_multi}}.  Default \code{NULL} (empirical mode).
#' @param constraint,cost_total,cost_fixed,cost_per_unit,cv_targets,n_min,method
#'   Design constraints and solver settings passed unchanged to
#'   \code{\link{stratify_multi}} for every candidate grid.
#' @param n_starts,max_iter Passed to \code{\link{stratify_multi}}.  Defaults
#'   \code{20} and \code{200}.
#' @param tol_gain Relative marginal-improvement threshold used to pick the
#'   recommended grid (default \code{0.02}, i.e. 2\%).
#' @param plot Logical.  Draw the \eqn{\Phi^*}-vs-\eqn{H_{\mathrm{total}}}
#'   curve with the recommended grid marked.  Default \code{TRUE}.
#' @param verbose Print progress per candidate.  Default \code{FALSE}.
#'
#' @return An object of class \code{"chooseH"}: a list with \code{table}
#'   (a data.frame with columns \code{H_label}, \code{H_total},
#'   \code{objective}, \code{optimality_gap}, \code{re_vs_coarsest} and
#'   \code{marginal_gain}), \code{recommended} (the chosen \code{H_vec}),
#'   \code{recommended_row}, and \code{results} (the list of fitted
#'   \code{"multistrata"} objects).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "independent")
#' ch  <- choose_H(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
#'                 n_total = 400, H_each = 2:4, n_starts = 5, plot = FALSE)
#' ch$table
#' ch$recommended
#' }
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{find_n}}
#' @export
choose_H <- function(data,
                     study_vars,
                     n_total,
                     H_list        = NULL,
                     H_each        = 2:5,
                     omega         = NULL,
                     distributions = NULL,
                     constraint    = c("n_total", "cost"),
                     cost_total    = NULL,
                     cost_fixed    = 0,
                     cost_per_unit = 1,
                     cv_targets    = NULL,
                     n_min         = 2L,
                     method        = c("cobyla", "slsqp"),
                     n_starts      = 20L,
                     max_iter      = 200L,
                     tol_gain      = 0.02,
                     plot          = TRUE,
                     verbose       = FALSE) {

  data       <- as.matrix(data)
  study_vars <- as.matrix(study_vars)
  constraint <- match.arg(constraint)
  method     <- match.arg(method)
  p <- ncol(data);  q <- ncol(study_vars)

  if (is.null(omega)) omega <- rep(1 / q, q)

  # Build candidate grid list
  if (is.null(H_list)) {
    H_list <- lapply(H_each, function(h) rep(as.integer(h), p))
  } else {
    if (!is.list(H_list) || !all(vapply(H_list, length, 1L) == p))
      stop("'H_list' must be a list of integer vectors, each of length p = ", p, ".")
    H_list <- lapply(H_list, as.integer)
  }

  # Order candidates by total number of strata (coarsest first)
  H_tot   <- vapply(H_list, prod, 1)
  ord     <- order(H_tot)
  H_list  <- H_list[ord]
  H_tot   <- H_tot[ord]
  n_cand  <- length(H_list)

  results  <- vector("list", n_cand)
  phi_vec  <- numeric(n_cand)
  gap_vec  <- numeric(n_cand)
  labels   <- character(n_cand)

  for (i in seq_len(n_cand)) {
    Hv <- H_list[[i]]
    labels[i] <- paste(Hv, collapse = "x")
    if (verbose)
      message(sprintf("  [choose_H] fitting H = %s (H_total = %d) ...",
                      labels[i], H_tot[i]))
    res <- stratify_multi(
      data = data, study_vars = study_vars, H_vec = Hv, n_total = n_total,
      omega = omega, constraint = constraint, cost_total = cost_total,
      cost_fixed = cost_fixed, cost_per_unit = cost_per_unit,
      cv_targets = cv_targets, n_min = n_min, method = method,
      n_starts = n_starts, max_iter = max_iter,
      distributions = distributions, verbose = FALSE)
    results[[i]] <- res
    phi_vec[i]   <- res$objective
    gap_vec[i]   <- res$optimality_gap
  }

  # Relative efficiency vs the coarsest grid, and marginal gain between
  # successive grids (fraction by which Phi* dropped from the previous one).
  re_vs_coarse  <- phi_vec[1L] / phi_vec
  marginal_gain <- c(NA_real_, -diff(phi_vec) / phi_vec[-n_cand])

  tab <- data.frame(
    H_label        = labels,
    H_total        = as.integer(H_tot),
    objective      = phi_vec,
    optimality_gap = gap_vec,
    re_vs_coarsest = re_vs_coarse,
    marginal_gain  = marginal_gain,
    stringsAsFactors = FALSE
  )

  # Elbow: first grid whose *next* step improves Phi* by less than tol_gain
  # (i.e. diminishing returns set in).  Fall back to the finest grid.
  rec_i <- n_cand
  for (i in seq_len(n_cand - 1L)) {
    if (marginal_gain[i + 1L] < tol_gain) { rec_i <- i; break }
  }

  out <- structure(
    list(
      table           = tab,
      recommended     = H_list[[rec_i]],
      recommended_row = rec_i,
      results         = results,
      tol_gain        = tol_gain
    ),
    class = "chooseH"
  )

  if (isTRUE(plot)) plot(out)
  out
}


#' Print method for chooseH
#' @param x A \code{"chooseH"} object.
#' @param ... Unused.
#' @export
print.chooseH <- function(x, ...) {
  cat("\n")
  cat("== choose_H(): number-of-strata trace ==========\n")
  disp <- x$table
  disp$objective      <- signif(disp$objective, 6)
  disp$optimality_gap <- sprintf("%.2f%%", 100 * disp$optimality_gap)
  disp$re_vs_coarsest <- sprintf("%.3f", disp$re_vs_coarsest)
  disp$marginal_gain  <- ifelse(is.na(disp$marginal_gain), "-",
                                sprintf("%.1f%%", 100 * disp$marginal_gain))
  print(disp, row.names = FALSE)
  cat("------------------------------------------------\n")
  cat(sprintf("  Recommended grid : %s  (H_total = %d)\n",
              paste(x$recommended, collapse = " x "), prod(x$recommended)))
  cat(sprintf("  Rule             : first grid whose next refinement gains\n"))
  cat(sprintf("                     < %.0f%% in Phi* (elbow).\n", 100 * x$tol_gain))
  cat("================================================\n\n")
  invisible(x)
}


#' Plot method for chooseH
#'
#' Draws the \eqn{\Phi^*}-vs-\eqn{H_{\mathrm{total}}} curve with the
#' recommended grid highlighted.
#'
#' @param x A \code{"chooseH"} object.
#' @param ... Further graphical parameters passed to \code{plot}.
#' @export
plot.chooseH <- function(x, ...) {
  tab <- x$table
  op  <- par(mar = c(4.6, 4.6, 3.2, 1.2))
  on.exit(par(op))

  plot(tab$H_total, tab$objective, type = "b", pch = 19, lwd = 2,
       col = "#4E79A7", bty = "l",
       xlab = expression("Total number of strata  " * H[total]),
       ylab = expression("Optimal composite variance  " * Phi^"*"),
       main = "choose_H(): precision vs. number of strata", ...)

  # Mark the recommended grid
  ri <- x$recommended_row
  points(tab$H_total[ri], tab$objective[ri], pch = 21, cex = 2.4,
         lwd = 2.5, col = "#E15759", bg = adjustcolor("#E15759", 0.25))
  text(tab$H_total[ri], tab$objective[ri],
       labels = sprintf("  recommended: %s", tab$H_label[ri]),
       pos = 4, col = "#E15759", font = 2, cex = 0.85, xpd = NA)

  # Annotate each point with its grid label
  text(tab$H_total, tab$objective, labels = tab$H_label,
       pos = 3, cex = 0.7, col = "grey35", xpd = NA)

  invisible(x)
}
