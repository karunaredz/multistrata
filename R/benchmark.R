# =============================================================================
# benchmark.R
#
# Fix vs original:
#   Called .round_to_total() directly, that function lives in neyman_alloc.R
#   and starts with a dot so it is an internal function.  Calling it from
#   benchmark.R works when sourced together but is fragile in a package.
#   Replaced with an inline Hamilton rounding step.
# =============================================================================


#' Benchmark Stratification Designs
#'
#' Compares the joint multivariate stratification result against proportional,
#' equal, and sequential-Neyman designs.  Operationalises Lemma 1 of
#' Reddy & Khan (forthcoming): sequential procedures are strictly sub-optimal.
#'
#' @param data         Matrix/data.frame of stratification variables (\eqn{N \times p}).
#' @param study_vars   Matrix/data.frame of study variables (\eqn{N \times q}).
#' @param H_vec        Integer vector: strata counts per variable.
#' @param n_total      Total sample size.
#' @param omega        Study-variable weight vector (length \eqn{q}).
#' @param result_joint A \code{"multistrata"} object.  If \code{NULL} the
#'   joint method is run internally.
#' @param n_starts     Multi-starts for the joint method if re-run.  Default 20.
#'
#' @return A \code{data.frame} with one row per method.
#'
#' @examples
#' set.seed(99)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "correlated")
#' bench <- benchmark_designs(
#'   data       = pop[, c("X1","X2")],
#'   study_vars = pop[, c("Y1","Y2")],
#'   H_vec      = c(3, 3),
#'   n_total    = 300,
#'   omega      = c(0.5, 0.5),
#'   n_starts   = 5
#' )
#' print(bench)
#'
#' @export
benchmark_designs <- function(data,
                               study_vars,
                               H_vec,
                               n_total,
                               omega        = NULL,
                               result_joint = NULL,
                               n_starts     = 20L) {

  data       <- as.matrix(data)
  study_vars <- as.matrix(study_vars)
  N  <- nrow(data);  p <- ncol(data);  q <- ncol(study_vars)

  if (is.null(omega)) omega <- rep(1 / q, q)
  n_min <- 2L

  # Run joint method if not supplied
  if (is.null(result_joint)) {
    message("Running joint stratify_multi ... (n_starts = ", n_starts, ")")
    result_joint <- stratify_multi(
      data = data, study_vars = study_vars, H_vec = H_vec,
      n_total = n_total, omega = omega, n_starts = n_starts, verbose = FALSE)
  }
  phi_joint <- result_joint$objective
  results   <- list()

  # (1) Joint
  results[["Joint (multistrata)"]] <- phi_joint

  # (2) Proportional allocation on joint boundaries
  ss_joint    <- result_joint$stratum_stats
  sigma2_cols <- grep("^sigma2_k", names(ss_joint), value = TRUE)
  sigma2_mat  <- as.matrix(ss_joint[, sigma2_cols])
  W_h         <- ss_joint$W_h
  N_h         <- ss_joint$N_h

  n_prop      <- pmax(n_min, round(W_h * n_total))
  # FIX: inline Hamilton rounding rather than calling internal .round_to_total
  n_prop      <- .bm_round_to_total(n_prop, n_total, n_min, N_h)
  phi_prop    <- composite_variance(sigma2_mat, W_h, n_prop, N_h, omega)
  results[["Proportional allocation"]] <- phi_prop

  # (3) Equal allocation
  H_total  <- prod(H_vec)
  n_equal  <- rep(floor(n_total / H_total), H_total)
  rem      <- n_total - sum(n_equal)
  if (rem > 0) n_equal[seq_len(rem)] <- n_equal[seq_len(rem)] + 1L
  n_equal  <- pmax(n_equal, n_min)
  phi_equal <- composite_variance(sigma2_mat, W_h, n_equal, N_h, omega)
  results[["Equal allocation"]] <- phi_equal

  # (4) Sequential: marginal quantile boundaries + Neyman allocation
  B_seq <- vector("list", p)
  for (j in seq_len(p)) {
    probs      <- seq(0, 1, length.out = H_vec[j] + 1L)[2:H_vec[j]]
    B_seq[[j]] <- quantile(data[, j], probs = probs, names = FALSE)
  }
  sv_seq       <- .compute_stratum_means_vars(data, study_vars, B_seq, H_vec)
  W_seq        <- sv_seq$counts / N
  N_seq        <- sv_seq$counts
  n_neyman_seq <- neyman_alloc(sv_seq$sigma2, W_seq, omega, n_total, N_seq, n_min)
  phi_seq      <- composite_variance(sv_seq$sigma2, W_seq, n_neyman_seq, N_seq, omega)
  results[["Sequential (quantile + Neyman)"]] <- phi_seq

  # Assemble
  phi_vec <- unlist(results)
  re_vec  <- phi_vec / phi_joint

  data.frame(
    method              = names(results),
    composite_variance  = round(phi_vec, 8),
    relative_efficiency = round(re_vec,  4),
    pct_gain_over_joint = round(100 * (re_vec - 1), 2),
    stringsAsFactors    = FALSE,
    row.names           = NULL
  )
}

# Internal Hamilton rounding for benchmark.R only.
# Mirrors .round_to_total() from neyman_alloc.R, includes the safety pass
# (Bug 5 fix) so sum(alloc) == target even when n_min is tight.
.bm_round_to_total <- function(raw, target, n_min, N_h) {
  H       <- length(raw)
  frac    <- raw - floor(raw)
  alloc   <- pmax(floor(raw), n_min)
  alloc   <- pmin(alloc, N_h)
  deficit <- target - sum(alloc)
  if (deficit > 0L) {
    ord <- order(frac, decreasing = TRUE)
    for (i in seq_len(min(as.integer(deficit), H))) {
      h <- ord[i]
      if (alloc[h] < N_h[h]) alloc[h] <- alloc[h] + 1L
    }
  } else if (deficit < 0L) {
    ord <- order(frac)
    rem <- 0L
    for (i in seq_len(H)) {
      if (rem >= abs(deficit)) break
      h <- ord[i]
      if (alloc[h] > n_min) { alloc[h] <- alloc[h] - 1L; rem <- rem + 1L }
    }
    # Safety pass: force-reduce from largest strata if still over-allocated
    remaining <- sum(alloc) - target
    if (remaining > 0L) {
      ord2 <- order(alloc, decreasing = TRUE)
      for (i in seq_len(H)) {
        if (remaining <= 0L) break
        h <- ord2[i]
        if (alloc[h] > n_min) { alloc[h] <- alloc[h] - 1L; remaining <- remaining - 1L }
      }
      if (remaining > 0L)
        warning(sprintf(
          ".bm_round_to_total: could not achieve sum = %d (n_min too tight); got %d.",
          target, as.integer(sum(alloc))))
    }
  }
  as.numeric(alloc)
}
