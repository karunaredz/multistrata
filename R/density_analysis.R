# =============================================================================
# density_analysis.R
#
# Downstream analysis helpers for density-only stratification.  These mirror the
# empirical helpers benchmark_designs(), find_n() and choose_H(), but work from
# analytic marginal densities rather than a realised population frame.
# =============================================================================


#' Benchmark Density-Based Stratification Designs
#'
#' Compares a joint density-based optimum against simple allocation and
#' boundaries-then-allocation alternatives.  This is the analytic, data-free
#' counterpart of \code{\link{benchmark_designs}} for results returned by
#' \code{\link{stratify_multi_density}}.
#'
#' @param distributions A named list of marginal density specifications.
#' @param H_vec Integer vector: strata counts per variable.
#' @param n_total Total sample size.  Defaults to the allocation total stored
#'   in \code{result_joint}, when supplied.
#' @param N Population size used for finite-population corrections.
#' @param ranges Optional support ranges, as in \code{\link{stratify_multi_density}}.
#' @param omega Study-variable weights.  Defaults to equal weights, or the
#'   weights stored in \code{result_joint}.
#' @param result_joint Optional \code{"multistrata"} density result.  If
#'   \code{NULL}, the joint method is run internally.
#' @param constraint,cost_total,cost_fixed,cost_per_unit,n_min Passed to
#'   \code{\link{stratify_multi_density}} / \code{\link{neyman_alloc}} where
#'   relevant.
#' @param n_starts,max_iter Passed to \code{\link{stratify_multi_density}} if
#'   the joint result is run internally.
#' @param subdivisions Maximum subdivisions passed to \code{stats::integrate}
#'   for density moment quadrature.
#'
#' @return A \code{data.frame} with one row per method.
#'
#' @examples
#' \donttest{
#' d <- list(X1 = list(family = "exp", rate = 1),
#'           X2 = list(family = "triang", a = 0, b = 2, c = 0))
#' res <- stratify_multi_density(d, H_vec = c(3, 3), n_total = 300,
#'                               N = 2000, ranges = list(c(0, 8), c(0, 2)),
#'                               n_starts = 3)
#' benchmark_density_designs(d, result_joint = res)
#' }
#'
#' @export
benchmark_density_designs <- function(distributions = NULL,
                                      H_vec         = NULL,
                                      n_total      = NULL,
                                      N            = NULL,
                                      ranges       = NULL,
                                      omega        = NULL,
                                      result_joint = NULL,
                                      constraint   = c("n_total", "cost"),
                                      cost_total   = NULL,
                                      cost_fixed   = 0,
                                      cost_per_unit = 1,
                                      n_min        = 2L,
                                      n_starts     = 20L,
                                      max_iter     = 100L,
                                      subdivisions = 200L) {

  constraint <- match.arg(constraint)
  meta <- .density_resolve_meta(distributions, H_vec, n_total, N, ranges,
                                omega, result_joint)
  distributions <- meta$distributions
  H_vec <- meta$H_vec
  n_total <- meta$n_total
  N <- meta$N
  ranges <- meta$ranges
  omega <- meta$omega

  if (is.null(result_joint)) {
    result_joint <- stratify_multi_density(
      distributions = distributions, H_vec = H_vec, n_total = n_total, N = N,
      ranges = ranges, omega = omega, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit, n_min = n_min, n_starts = n_starts,
      max_iter = max_iter, subdivisions = subdivisions, verbose = FALSE)
  } else if (!identical(result_joint$.meta$mode, "density")) {
    stop("'result_joint' must be a density result from stratify_multi_density().")
  }

  phi_joint <- result_joint$objective
  ss <- result_joint$stratum_stats
  mats <- .density_extract_moments(ss, length(omega))
  W_h <- ss$W_h
  N_h <- ss$N_h

  results <- list("Joint (multistrata density)" = phi_joint)

  # Proportional allocation on joint boundaries.
  n_prop <- .bm_round_to_total(W_h * n_total, as.integer(round(n_total)),
                               n_min, N_h)
  results[["Proportional allocation"]] <-
    composite_variance(mats$sigma2, W_h, n_prop, N_h, omega)

  # Equal allocation on joint boundaries.
  n_equal <- .bm_round_to_total(rep(n_total / length(W_h), length(W_h)),
                                as.integer(round(n_total)), n_min, N_h)
  results[["Equal allocation"]] <-
    composite_variance(mats$sigma2, W_h, n_equal, N_h, omega)

  # Sequential analogue: marginal density quantile boundaries + Neyman.
  B_seq <- .density_quantile_boundaries(distributions, H_vec, ranges)
  asm_seq <- .density_assemble(distributions, B_seq, H_vec, ranges, N,
                               subdivisions)
  n_seq <- neyman_alloc(asm_seq$sigma2, asm_seq$W_h, omega,
                        n_total = n_total, N_h = asm_seq$N_h, n_min = n_min,
                        constraint = constraint, cost_total = cost_total,
                        cost_fixed = cost_fixed,
                        cost_per_unit = cost_per_unit, mu = asm_seq$mu)
  results[["Sequential (density quantile + Neyman)"]] <-
    composite_variance(asm_seq$sigma2, asm_seq$W_h, n_seq, asm_seq$N_h, omega)

  phi_vec <- unlist(results)
  re_vec <- phi_vec / phi_joint
  data.frame(
    method = names(results),
    composite_variance = round(phi_vec, 8),
    relative_efficiency = round(re_vec, 4),
    pct_gain_over_joint = round(100 * (re_vec - 1), 2),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Smallest Sample Size for a Density-Based Design
#'
#' Finds the smallest total sample size that satisfies coefficient-of-variation
#' targets for a density-only stratification.  With \code{reoptimize = FALSE}
#' the current density boundaries are reused; with \code{TRUE}, boundaries are
#' re-optimised at each bisection probe.
#'
#' @inheritParams benchmark_density_designs
#' @param cv_targets Numeric vector of maximum acceptable CV values.
#' @param result Optional density \code{"multistrata"} object supplying the
#'   current boundaries.
#' @param n_min Minimum sample size per stratum.
#' @param n_range Optional integer bracket \code{c(n_lo, n_hi)}.  Defaults to
#'   \code{c(H_total * n_min, N)}.
#' @param reoptimize Logical; re-optimise density boundaries at every trial
#'   sample size.  Default \code{FALSE}.
#' @param n_starts,max_iter Passed to \code{\link{stratify_multi_density}} when
#'   \code{reoptimize = TRUE}.
#' @param subdivisions Maximum subdivisions passed to \code{stats::integrate}
#'   for density moment quadrature.
#' @param verbose Print bisection probes.
#'
#' @return A \code{"findN"} object.
#' @export
find_n_density <- function(distributions = NULL,
                           cv_targets,
                           H_vec         = NULL,
                           N             = NULL,
                           ranges        = NULL,
                           omega         = NULL,
                           result        = NULL,
                           n_min         = 2L,
                           n_range       = NULL,
                           reoptimize    = FALSE,
                           n_starts      = 20L,
                           max_iter      = 100L,
                           subdivisions  = 200L,
                           verbose       = FALSE) {

  meta <- .density_resolve_meta(distributions, H_vec, NULL, N, ranges, omega,
                                result)
  distributions <- meta$distributions
  H_vec <- meta$H_vec
  N <- meta$N
  ranges <- meta$ranges
  omega <- meta$omega
  q <- length(omega)

  if (missing(cv_targets) || is.null(cv_targets))
    stop("'cv_targets' must be supplied.")
  if (length(cv_targets) != q)
    stop("'cv_targets' must have length ", q, ".")
  if (any(cv_targets <= 0))
    stop("All elements of 'cv_targets' must be strictly positive.")

  H_total <- prod(H_vec)
  if (is.null(n_range)) n_range <- c(H_total * n_min, N)
  n_range <- as.integer(round(n_range))
  n_lo <- max(n_range[1L], H_total * n_min)
  n_hi <- min(n_range[2L], N)
  if (n_hi <= n_lo)
    stop("Empty search bracket: n_range = (", n_lo, ", ", n_hi, ").")

  B <- if (!is.null(result)) result$boundaries else
    .density_quantile_boundaries(distributions, H_vec, ranges)

  if (!reoptimize) {
    asm <- .density_assemble(distributions, B, H_vec, ranges, N, subdivisions)
    feasible_at <- function(n) {
      alloc <- suppressWarnings(neyman_alloc(
        sigma2 = asm$sigma2, W_h = asm$W_h, omega = omega, n_total = n,
        N_h = asm$N_h, n_min = n_min, cv_targets = cv_targets, mu = asm$mu))
      cv <- attr(alloc, "cv_achieved")
      ok <- isTRUE(attr(alloc, "cv_feasible")) &&
        all(cv <= cv_targets + 5e-4)
      list(feasible = ok, alloc = alloc, cv = cv)
    }
  } else {
    feasible_at <- function(n) {
      res <- suppressWarnings(stratify_multi_density(
        distributions = distributions, H_vec = H_vec, n_total = n, N = N,
        ranges = ranges, omega = omega, cv_targets = cv_targets,
        n_min = n_min, n_starts = n_starts, max_iter = max_iter,
        subdivisions = subdivisions, verbose = FALSE))
      cv <- res$cv_achieved
      ok <- !is.null(cv) && all(cv <= cv_targets + 5e-4)
      list(feasible = ok, alloc = res$allocation, cv = cv, result = res)
    }
  }

  trace <- data.frame(n = integer(0), feasible = logical(0))
  probe <- function(n) {
    r <- feasible_at(n)
    trace[nrow(trace) + 1L, ] <<- list(as.integer(n), r$feasible)
    if (verbose)
      message(sprintf("  [find_n_density] n = %6d -> %s",
                      as.integer(n), if (r$feasible) "feasible" else "infeasible"))
    r
  }

  r_hi <- probe(n_hi)
  if (!r_hi$feasible) {
    warning("CV targets are unattainable within n_range.")
    return(.make_findN(NA_integer_, r_hi$alloc, r_hi$cv, cv_targets, FALSE,
                       omega = omega, B = B, H_vec = H_vec,
                       n_range = c(n_lo, n_hi),
                       trace = trace[order(trace$n), , drop = FALSE]))
  }

  r_lo <- probe(n_lo)
  if (r_lo$feasible) {
    best_n <- n_lo
    best_r <- r_lo
  } else {
    lo <- n_lo
    hi <- n_hi
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      rm <- probe(mid)
      if (rm$feasible) hi <- mid else lo <- mid
    }
    best_n <- hi
    best_r <- probe(hi)
  }

  if (!reoptimize) {
    phi <- composite_variance(asm$sigma2, asm$W_h, best_r$alloc,
                              asm$N_h, omega)
    B_out <- B
  } else {
    phi <- best_r$result$objective
    B_out <- best_r$result$boundaries
  }

  .make_findN(as.integer(best_n), best_r$alloc, best_r$cv, cv_targets, TRUE,
              omega = omega, B = B_out, H_vec = H_vec, phi = phi,
              n_range = c(n_lo, n_hi),
              trace = trace[order(trace$n), , drop = FALSE])
}


#' Choose the Number of Strata for Density-Based Stratification
#'
#' Analytic density-mode counterpart of \code{\link{choose_H}}.  It traces the
#' optimum composite variance over candidate strata grids and recommends the
#' first elbow where the next refinement gives less than \code{tol_gain}
#' improvement.
#'
#' @inheritParams benchmark_density_designs
#' @param H_list Optional list of integer \code{H_vec} candidates.
#' @param H_each Used to build \code{H_list} as \code{rep(h, p)} when
#'   \code{H_list} is \code{NULL}.
#' @param cv_targets Optional CV targets passed to
#'   \code{\link{stratify_multi_density}}.
#' @param subdivisions Maximum subdivisions passed to \code{stats::integrate}
#'   for density moment quadrature.
#' @param tol_gain Elbow threshold.
#' @param plot Logical; draw the trace.
#' @param verbose Logical; print progress.
#'
#' @return A \code{"chooseH"} object.
#' @export
choose_H_density <- function(distributions,
                             n_total,
                             N,
                             H_list        = NULL,
                             H_each        = 2:5,
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
                             subdivisions  = 200L,
                             tol_gain      = 0.02,
                             plot          = TRUE,
                             verbose       = FALSE) {

  constraint <- match.arg(constraint)
  p <- length(distributions)
  if (is.null(omega)) omega <- rep(1 / p, p)
  if (length(omega) != p) stop("'omega' must have length ", p, ".")

  if (is.null(ranges)) ranges <- .density_default_ranges(distributions)
  if (length(ranges) != p) stop("'ranges' must have length ", p, ".")

  if (is.null(H_list)) {
    H_list <- lapply(H_each, function(h) rep(as.integer(h), p))
  } else {
    if (!is.list(H_list) || !all(vapply(H_list, length, 1L) == p))
      stop("'H_list' must be a list of integer vectors, each of length p = ", p, ".")
    H_list <- lapply(H_list, as.integer)
  }

  H_tot <- vapply(H_list, prod, 1)
  ord <- order(H_tot)
  H_list <- H_list[ord]
  H_tot <- H_tot[ord]
  n_cand <- length(H_list)

  results <- vector("list", n_cand)
  phi_vec <- numeric(n_cand)
  gap_vec <- numeric(n_cand)
  labels <- character(n_cand)

  for (i in seq_len(n_cand)) {
    Hv <- H_list[[i]]
    labels[i] <- paste(Hv, collapse = "x")
    if (verbose)
      message(sprintf("  [choose_H_density] fitting H = %s ...", labels[i]))
    res <- stratify_multi_density(
      distributions = distributions, H_vec = Hv, n_total = n_total, N = N,
      ranges = ranges, omega = omega, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit, cv_targets = cv_targets, n_min = n_min,
      n_starts = n_starts, max_iter = max_iter, subdivisions = subdivisions,
      verbose = FALSE)
    results[[i]] <- res
    phi_vec[i] <- res$objective
    gap_vec[i] <- res$optimality_gap
  }

  re_vs_coarse <- phi_vec[1L] / phi_vec
  marginal_gain <- c(NA_real_, -diff(phi_vec) / phi_vec[-n_cand])
  tab <- data.frame(
    H_label = labels,
    H_total = as.integer(H_tot),
    objective = phi_vec,
    optimality_gap = gap_vec,
    re_vs_coarsest = re_vs_coarse,
    marginal_gain = marginal_gain,
    stringsAsFactors = FALSE
  )

  rec_i <- n_cand
  for (i in seq_len(n_cand - 1L)) {
    if (marginal_gain[i + 1L] < tol_gain) {
      rec_i <- i
      break
    }
  }

  out <- structure(list(table = tab, recommended = H_list[[rec_i]],
                        recommended_row = rec_i, results = results,
                        tol_gain = tol_gain),
                   class = "chooseH")
  if (isTRUE(plot)) plot(out)
  out
}


#' Omega Sweep for Density-Based Stratification
#'
#' Re-optimises a density-based design over a grid of priority weights for two
#' variables, returning objective and CV diagnostics for each grid point.
#'
#' @inheritParams benchmark_density_designs
#' @param pair Integer vector of length 2 giving the variables whose relative
#'   weights are swept. Defaults to the first two variables.
#' @param grid Numeric vector of shares assigned to \code{pair[1]} within the
#'   pair's total weight.
#' @param cv_targets Optional CV targets passed to
#'   \code{\link{stratify_multi_density}}.
#' @param subdivisions Maximum subdivisions passed to \code{stats::integrate}
#'   for density moment quadrature.
#'
#' @return A \code{data.frame} with one row per grid point and columns
#'   \code{w}, \code{phi}, \code{cv1}, \code{cv2} and \code{gap}, matching
#'   \code{\link{omega_sweep}}. The swept indices and names are attached as the
#'   \code{"pair"} and \code{"pair_names"} attributes.
#'
#' @seealso \code{\link{omega_sweep}} for the same sweep on observed data.
#' @export
omega_sweep_density <- function(distributions,
                                H_vec,
                                n_total,
                                N,
                                ranges        = NULL,
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
                                subdivisions  = 200L) {

  constraint <- match.arg(constraint)
  p <- length(distributions)
  if (is.null(omega)) omega <- rep(1 / p, p)
  if (length(pair) != 2L || any(pair < 1L) || any(pair > p) || pair[1] == pair[2])
    stop("'pair' must contain two different variable indices.")
  if (length(omega) != p) stop("'omega' must have length ", p, ".")
  if (is.null(ranges)) ranges <- .density_default_ranges(distributions)

  total_pair <- sum(omega[pair])
  out <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    om <- omega
    om[pair[1]] <- grid[i] * total_pair
    om[pair[2]] <- (1 - grid[i]) * total_pair
    om <- om / sum(om)
    r <- stratify_multi_density(
      distributions = distributions, H_vec = H_vec, n_total = n_total, N = N,
      ranges = ranges, omega = om, constraint = constraint,
      cost_total = cost_total, cost_fixed = cost_fixed,
      cost_per_unit = cost_per_unit, cv_targets = cv_targets, n_min = n_min,
      n_starts = n_starts, max_iter = max_iter, subdivisions = subdivisions,
      verbose = FALSE)
    mats <- .density_extract_moments(r$stratum_stats, p)
    cv <- as.numeric(design_cv(mats$sigma2, mats$mu, r$stratum_stats$W_h,
                               r$allocation, r$stratum_stats$N_h))
    out[[i]] <- data.frame(
      w = om[pair[1]],
      phi = r$objective,
      cv1 = 100 * cv[pair[1]],
      cv2 = 100 * cv[pair[2]],
      gap = 100 * r$optimality_gap,
      stringsAsFactors = FALSE
    )
  }
  tab <- do.call(rbind, out)
  d_names <- names(distributions)
  if (is.null(d_names)) d_names <- paste0("X", seq_len(p))
  attr(tab, "pair")       <- as.integer(pair)
  attr(tab, "pair_names") <- d_names[pair]
  tab
}


# Internal density-analysis helpers ------------------------------------------

.density_resolve_meta <- function(distributions, H_vec, n_total, N, ranges,
                                  omega, result = NULL) {
  if (!is.null(result)) {
    if (!inherits(result, "multistrata") ||
        !identical(result$.meta$mode, "density"))
      stop("'result' must be a density 'multistrata' object.")
    m <- result$.meta
    if (missing(distributions) || is.null(distributions))
      distributions <- m$distributions
    if (is.null(H_vec)) H_vec <- m$H_vec
    if (is.null(n_total)) n_total <- sum(result$allocation)
    if (is.null(N)) N <- m$N
    if (is.null(ranges)) ranges <- result$ranges
    if (is.null(omega)) omega <- m$omega
  }
  if (!is.list(distributions))
    stop("'distributions' must be supplied as a list.")
  p <- length(distributions)
  if (is.null(H_vec)) stop("'H_vec' must be supplied.")
  if (length(H_vec) != p) stop("'H_vec' must have length ", p, ".")
  if (is.null(N)) stop("'N' must be supplied.")
  if (is.null(ranges)) ranges <- .density_default_ranges(distributions)
  if (length(ranges) != p) stop("'ranges' must have length ", p, ".")
  if (is.null(omega)) omega <- rep(1 / p, p)
  if (length(omega) != p) stop("'omega' must have length ", p, ".")
  if (abs(sum(omega) - 1) > 1e-8) stop("'omega' must sum to 1.")
  if (!is.null(n_total)) n_total <- as.integer(round(n_total))
  list(distributions = distributions, H_vec = as.integer(H_vec),
       n_total = n_total, N = as.integer(round(N)), ranges = ranges,
       omega = omega)
}

.density_default_ranges <- function(distributions) {
  lapply(seq_along(distributions), function(j) {
    qf <- .build_quantile(distributions[[j]], 0, 1)
    c(qf(1e-6), qf(1 - 1e-6))
  })
}

.density_quantile_boundaries <- function(distributions, H_vec, ranges) {
  lapply(seq_along(H_vec), function(j) {
    qf <- .build_quantile(distributions[[j]], ranges[[j]][1L], ranges[[j]][2L])
    probs <- seq(0, 1, length.out = H_vec[j] + 1L)[2:H_vec[j]]
    as.numeric(qf(probs))
  })
}

.density_assemble <- function(distributions, B, H_vec, ranges, N,
                              subdivisions = 200L) {
  margins <- .density_margins(distributions, B, H_vec, ranges, subdivisions)
  .assemble_from_margins(margins, H_vec, N)
}

.density_extract_moments <- function(ss, q) {
  # Same column layout as the data mode; shared implementation lives in
  # omega_sweep.R so both branches cannot drift apart.
  .extract_stratum_moments(ss, q)
}
