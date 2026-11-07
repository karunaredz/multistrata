# =============================================================================
# survey_tools.R
#
# Survey-ecosystem plumbing: score an arbitrary frame with the optimal strata
# (assign_strata) and export a ready-to-analyse design object for the survey
# package / srvyr (as_svydesign).  These remove the "last-mile" friction of
# taking a multistrata design into production estimation.
#
# Reference:
#   Reddy, K. G. & Khan, M. G. M. Joint Optimum Stratification and Sample
#   Allocation for Multivariate Surveys (upcoming methods paper).
# =============================================================================


#' Assign Strata to a New Frame
#'
#' Applies the optimal boundaries from a fitted \code{"multistrata"} object
#' to score any set of units (the estimation frame, a fresh census extract, or
#' hold-out data) into their strata.  This is the operational counterpart of
#' \code{\link{stratify_multi}}: the design is estimated once, then reused to
#' label every unit that must be sampled or weighted.
#'
#' @param result A \code{"multistrata"} object from
#'   \code{\link{stratify_multi}}.
#' @param newdata A \code{data.frame} or \code{matrix} of stratification
#'   variables for the units to be scored.  If it carries the same column
#'   names as the variables used to fit \code{result}, those columns are
#'   selected (in the correct order) automatically; otherwise it must have
#'   exactly \eqn{p} columns, taken in the original order.
#' @param labels Logical.  Return a \code{factor} of stratum labels
#'   (\code{TRUE}, default) or an integer stratum index (\code{FALSE}).
#'
#' @return A vector of length \code{nrow(newdata)}: a factor of stratum labels
#'   (\code{labels = TRUE}) or integer stratum indices (\code{labels = FALSE}).
#'   The attribute \code{"index_matrix"} holds the per-variable marginal
#'   stratum indices (\eqn{n \times p}, one-based).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "independent")
#' res <- stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
#'                       H_vec = c(3, 3), n_total = 300, n_starts = 5)
#' s <- assign_strata(res, pop[, c("X1", "X2")])
#' table(s)
#' }
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{as_svydesign}}
#' @export
assign_strata <- function(result, newdata, labels = TRUE) {

  if (!inherits(result, "multistrata"))
    stop("'result' must be a 'multistrata' object from stratify_multi().")

  m      <- result$.meta
  p      <- m$p
  H_vec  <- m$H_vec
  x_nm   <- m$x_names

  nd <- .match_strat_columns(newdata, x_nm, p)

  strat_id <- .assign_strata(nd, result$boundaries, H_vec)

  # Per-variable marginal indices (one-based) for optional inspection
  idx_mat <- matrix(0L, nrow = nrow(nd), ncol = p,
                    dimnames = list(NULL, x_nm))
  for (j in seq_len(p)) {
    breaks       <- c(-Inf, result$boundaries[[j]], Inf)
    ij           <- as.integer(cut(nd[, j], breaks = breaks,
                                   right = FALSE, labels = FALSE))
    idx_mat[, j] <- pmax(1L, pmin(H_vec[j], ij))
  }

  out <- if (labels) {
    lv <- .stratum_index_labels(H_vec)
    factor(lv[strat_id], levels = lv)
  } else {
    strat_id
  }
  attr(out, "index_matrix") <- idx_mat
  out
}


#' Export a multistrata Design to the survey Package
#'
#' Builds a \code{survey::svydesign} object (or, with \code{frame_only = TRUE},
#' the augmented data frame) from a fitted \code{"multistrata"} object and a
#' set of sampled units, ready for design-based estimation with the
#' \pkg{survey} package or \pkg{srvyr}.  Each sampled unit is scored into its
#' stratum, assigned a design (inverse-inclusion-probability) weight
#' \eqn{w_i = N_{h} / n_{h}}, and given a finite-population correction
#' \eqn{N_{h}}.
#'
#' @param result A \code{"multistrata"} object from
#'   \code{\link{stratify_multi}}.
#' @param sample_data A \code{data.frame} of the sampled units.  Must contain
#'   the stratification variables (matched by name to those used to fit
#'   \code{result}, or as the first \eqn{p} columns) and typically the
#'   observed study variables as well.
#' @param use_planned Logical.  Use the planned per-stratum allocation
#'   \code{result$allocation} as \eqn{n_h} in the design weights
#'   (\code{TRUE}), or the realised per-stratum counts in \code{sample_data}
#'   (\code{FALSE}, default).  The realised counts give exact design weights
#'   for the sample actually drawn.
#' @param ids Cluster/PSU formula passed to \code{survey::svydesign}.  Default
#'   \code{~1} (element sampling, no clustering).
#' @param frame_only Logical.  Return only the augmented data frame (with the
#'   \code{.stratum}, \code{.N_h}, \code{.n_h} and \code{.weight} columns)
#'   instead of a \code{svydesign} object.  Useful for \pkg{srvyr} pipelines
#'   or when \pkg{survey} is not installed.  Default \code{FALSE}.
#'
#' @return A \code{survey.design} object, or (when \code{frame_only = TRUE})
#'   the augmented \code{data.frame}.
#'
#' @details
#' The added columns are: \code{.stratum} (factor of stratum labels),
#' \code{.N_h} (stratum population size, the finite-population correction),
#' \code{.n_h} (stratum sample size), and \code{.weight} (design weight
#' \eqn{N_h/n_h}).  With \pkg{srvyr} the frame can be turned into a survey
#' object via \code{srvyr::as_survey_design(frame, strata = .stratum,
#' weights = .weight, fpc = .N_h)}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' pop <- sim_population(N = 4000, p = 2, q = 2, scenario = "independent")
#' res <- stratify_multi(pop[, c("X1", "X2")], pop[, c("Y1", "Y2")],
#'                       H_vec = c(3, 3), n_total = 300, n_starts = 5)
#' # draw a stratified sample using the optimal allocation
#' s   <- assign_strata(res, pop[, c("X1", "X2")], labels = FALSE)
#' take <- unlist(lapply(seq_along(res$allocation), function(h) {
#'   idx <- which(s == h)
#'   sample(idx, min(length(idx), res$allocation[h]))
#' }))
#' frame <- as_svydesign(res, pop[take, ], frame_only = TRUE)
#' head(frame)
#' }
#'
#' @seealso \code{\link{assign_strata}}, \code{\link{stratify_multi}}
#' @export
as_svydesign <- function(result,
                         sample_data,
                         use_planned = FALSE,
                         ids         = ~1,
                         frame_only  = FALSE) {

  if (!inherits(result, "multistrata"))
    stop("'result' must be a 'multistrata' object from stratify_multi().")

  sample_data <- as.data.frame(sample_data)
  m     <- result$.meta
  H_vec <- m$H_vec
  ss    <- result$stratum_stats
  labs  <- ss$stratum

  # Score each sampled unit into its stratum (label + integer id)
  sid   <- as.integer(assign_strata(result, sample_data, labels = FALSE))
  strat <- factor(labs[sid], levels = labs)

  # Stratum population sizes (N_h) and sample sizes (n_h)
  N_h_by <- setNames(ss$N_h, labs)
  if (use_planned) {
    n_h_by <- setNames(as.integer(round(result$allocation)), labs)
  } else {
    n_h_by <- table(factor(strat, levels = labs))
    n_h_by <- setNames(as.integer(n_h_by), labs)
  }

  N_h_i <- as.numeric(N_h_by[as.character(strat)])
  n_h_i <- as.numeric(n_h_by[as.character(strat)])
  n_h_i[n_h_i < 1] <- NA_real_           # guard against divide-by-zero
  w_i   <- N_h_i / n_h_i

  frame <- sample_data
  frame$.stratum <- strat
  frame$.N_h     <- N_h_i
  frame$.n_h     <- n_h_i
  frame$.weight  <- w_i

  if (frame_only) return(frame)

  if (!requireNamespace("survey", quietly = TRUE))
    stop("The 'survey' package is required for as_svydesign(). ",
         "Install it, or call with frame_only = TRUE to obtain the ",
         "augmented data frame for use with srvyr or other tools.")

  survey::svydesign(
    ids     = ids,
    strata  = ~.stratum,
    weights = ~.weight,
    fpc     = ~.N_h,
    data    = frame
  )
}


# .match_strat_columns --------------------------------------------------------
# Select and order the p stratification columns from user data, matching by
# name when possible and otherwise falling back to positional use.
.match_strat_columns <- function(newdata, x_nm, p) {
  if (is.null(dim(newdata)))
    newdata <- matrix(newdata, ncol = 1L)
  cn <- colnames(newdata)
  if (!is.null(cn) && all(x_nm %in% cn)) {
    nd <- as.matrix(as.data.frame(newdata)[, x_nm, drop = FALSE])
  } else {
    if (ncol(newdata) != p)
      stop("'newdata' must have ", p, " columns (the stratification ",
           "variables ", paste(x_nm, collapse = ", "),
           "), or carry those column names.")
    nd <- as.matrix(newdata)
  }
  storage.mode(nd) <- "double"
  nd
}
