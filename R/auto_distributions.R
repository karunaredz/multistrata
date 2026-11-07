# =============================================================================
# auto_distributions.R
#
# Automated marginal family selection for the parametric mode
# ("auto-parametric" mode).  For each stratification variable X_j, candidate
# families are fitted by maximum likelihood (or closed-form / moment
# estimators where standard) and the family minimising AIC is selected.
# The winning specification is returned in exactly the list format accepted
# by stratify_multi(distributions = ...), so the analytic CDF/PDF/quantile
# machinery of stratum_integrals.R applies unchanged.
#
# Two engines:
#   "builtin", self-contained: MASS::fitdistr for norm/lnorm/gamma/weibull/
#                logis/cauchy (closed forms for exp/unif/pareto, moments for
#                the rescaled beta), log-likelihoods evaluated with the
#                package's own .build_pdf() so that the selected model is
#                *identical* to what the optimiser will use.
#   "fitverse", delegates the shared-parameterisation families (norm, lnorm,
#                gamma, exp, weibull, cauchy, logis, unif) to the FitVerse
#                package (Reddy & Khan, 2026), the data-characterisation
#                layer of the OptiStrata ecosystem, when it is installed.
#                FitVerse's Pareto/beta/Student-t use different
#                parameterisations (actuar/location-scale), so those
#                candidates are always handled by the builtin engine.
#   "auto", builtin candidates always; if FitVerse is installed its
#                best fit competes on AIC computed on the same data.
# =============================================================================


#' Automated Marginal Distribution Selection for the Parametric Mode
#'
#' Fits a set of candidate parametric families to each stratification
#' variable and returns the AIC-best specification per variable, in the
#' \code{distributions} list format accepted by \code{\link{stratify_multi}}.
#' Called internally when \code{stratify_multi(distributions = "auto")}.
#'
#' @param data Matrix or data.frame of stratification variables (columns).
#' @param families Character vector of candidate families.  Default:
#'   \code{c("lnorm", "gamma", "weibull", "norm", "exp", "logis", "cauchy",
#'   "pareto", "unif", "beta")}.  Families with positive support are skipped
#'   automatically for variables containing non-positive values.
#' @param engine One of \code{"auto"} (default), \code{"builtin"},
#'   \code{"fitverse"}.  See Details.
#' @param max_fit_n Variables longer than this are subsampled (with a fixed
#'   seed) for fitting speed.  Default \code{8000}.
#' @param verbose Print the selected family per variable.  Default
#'   \code{FALSE}.
#'
#' @details
#' The builtin engine evaluates every candidate's log-likelihood with the
#' package's own density constructors, so model selection is performed on
#' exactly the object the optimiser later uses.  When FitVerse is installed
#' (\code{engine = "auto"} or \code{"fitverse"}), the families whose
#' parameterisations coincide with base R (norm, lnorm, gamma, exp, weibull,
#' cauchy, logis, unif) are additionally fitted through
#' \code{FitVerse::fitverse()} and the overall AIC winner is retained.
#'
#' @return A named list of length \code{ncol(data)} suitable for
#'   \code{stratify_multi(distributions = )}.  Attribute \code{"fit_table"}
#'   holds a data.frame with one row per (variable, family) giving the
#'   estimated parameters, log-likelihood, AIC, and the engine used;
#'   attribute \code{"chosen"} gives the winning family per variable.
#'
#' @examples
#' set.seed(1)
#' pop <- sim_population(N = 3000, p = 2, q = 2, scenario = "skewed")
#' dl  <- fit_distributions(pop[, c("X1", "X2")], verbose = TRUE)
#' attr(dl, "chosen")
#' \donttest{
#' res <- stratify_multi(pop[, c("X1","X2")], pop[, c("Y1","Y2")],
#'                       H_vec = c(3, 3), n_total = 300, n_starts = 5,
#'                       distributions = "auto")
#' }
#'
#' @references
#' Reddy, K. G. & Khan, M. G. M. (2026). FitVerse: Comprehensive Distribution
#' Fitting, Parameter Estimation and Visualisation. R package.
#'
#' @export
fit_distributions <- function(data,
                              families = c("lnorm", "gamma", "weibull",
                                           "norm", "exp", "logis", "cauchy",
                                           "pareto", "unif", "beta"),
                              engine   = c("auto", "builtin", "fitverse"),
                              max_fit_n = 8000L,
                              verbose  = FALSE) {

  engine <- match.arg(engine)
  data   <- as.matrix(data)
  p      <- ncol(data)
  x_names <- if (!is.null(colnames(data))) colnames(data)
             else paste0("X", seq_len(p))

  fv_ok <- engine %in% c("auto", "fitverse") &&
           requireNamespace("FitVerse", quietly = TRUE)
  if (engine == "fitverse" && !fv_ok)
    stop("engine = 'fitverse' requires the FitVerse package. ",
         "Install it or use engine = 'builtin'.")

  out   <- vector("list", p)
  rows  <- list()
  chosen <- character(p)

  # --- RNG hygiene -----------------------------------------------------------
  # Columns longer than 'max_fit_n' are subsampled below under a fixed seed, so
  # that the AIC family selection is stable from run to run instead of
  # wobbling with whatever seed the caller happened to set.
  #
  # set.seed() writes to the global stream, however, and doing that without
  # putting it back had two unwanted effects. It changed the results of
  # whatever the caller did next; and because
  # stratify_multi(distributions = "auto") calls this function *before* it
  # draws its multi-start initial boundaries, it silently overrode the user's
  # own set.seed() for the whole optimisation that followed -- two different
  # user seeds produced byte-identical designs, so the design looked
  # reproducible while not actually being under the user's control.
  #
  # The stream is therefore saved here and restored on exit, and the RNG is
  # only touched at all when some column really is subsampled. To make the
  # subsample itself reproducible from the caller's seed instead, pass
  # max_fit_n >= nrow(data) so that no subsampling happens.
  if (nrow(data) > max_fit_n) {
    .seed_existed <- exists(".Random.seed", envir = globalenv(),
                            inherits = FALSE)
    .old_seed <- if (.seed_existed)
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
      if (.seed_existed)
        assign(".Random.seed", .old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
  }

  for (j in seq_len(p)) {
    x <- data[, j]
    x <- x[is.finite(x)]
    if (length(x) > max_fit_n) {
      set.seed(20260612L + j)
      x <- sample(x, max_fit_n)
    }
    a_j <- min(x); b_j <- max(x)

    cands <- .fit_candidates_builtin(x, families, a_j, b_j)
    for (cd in cands) rows[[length(rows) + 1L]] <-
      data.frame(variable = x_names[j], family = cd$family,
                 params = .fmt_params(cd$spec), logLik = cd$ll,
                 AIC = cd$aic, engine = "builtin",
                 stringsAsFactors = FALSE)

    if (fv_ok) {
      fv <- .fit_via_fitverse(x, families)
      if (!is.null(fv)) {
        rows[[length(rows) + 1L]] <-
          data.frame(variable = x_names[j], family = fv$family,
                     params = .fmt_params(fv$spec), logLik = fv$ll,
                     AIC = fv$aic, engine = "FitVerse",
                     stringsAsFactors = FALSE)
        cands <- c(cands, list(fv))
      }
    }

    if (!length(cands))
      stop("No candidate family could be fitted for variable '",
           x_names[j], "'. Supply 'distributions' manually.")

    best <- cands[[which.min(vapply(cands, function(cd) cd$aic, 1))]]
    out[[j]]  <- best$spec
    chosen[j] <- best$family
    if (verbose)
      message(sprintf("  [auto] %s ~ %s (AIC = %.1f, %s)",
                      x_names[j], best$family, best$aic, best$engine))
  }

  names(out) <- x_names
  ft <- do.call(rbind, rows)
  ft <- ft[order(ft$variable, ft$AIC), ]
  attr(out, "fit_table") <- ft
  attr(out, "chosen")    <- stats::setNames(chosen, x_names)
  out
}


# .fit_candidates_builtin -----------------------------------------------------
# Fit each requested family with a standard estimator; score with the
# package's own .build_pdf so selection matches the optimiser's model.
.fit_candidates_builtin <- function(x, families, a_j, b_j) {
  n   <- length(x)
  pos <- all(x > 0)
  out <- list()

  add <- function(family, spec, k) {
    pdf <- tryCatch(.build_pdf(spec, a_j, b_j), error = function(e) NULL)
    if (is.null(pdf)) return()
    dv <- tryCatch(pdf(x), error = function(e) NULL)
    if (is.null(dv) || any(!is.finite(dv)) || any(dv < 0)) return()
    if (mean(dv <= 0) > 0.001) return()               # support mismatch
    ll  <- sum(log(pmax(dv, 1e-300)))
    out[[length(out) + 1L]] <<- list(family = family, spec = spec,
                                     ll = ll, aic = 2 * k - 2 * ll,
                                     engine = "builtin")
  }

  fd <- function(dens) tryCatch(
    suppressWarnings(MASS::fitdistr(x, dens)), error = function(e) NULL)

  for (fam in families) {
    switch(fam,
      "norm" = {
        m <- mean(x); s <- sqrt(mean((x - m)^2))
        add("norm", list(family = "norm", mean = m, sd = s), 2L)
      },
      "lnorm" = if (pos) {
        ml <- mean(log(x)); sl <- sqrt(mean((log(x) - ml)^2))
        add("lnorm", list(family = "lnorm", meanlog = ml, sdlog = sl), 2L)
      },
      "exp" = if (pos)
        add("exp", list(family = "exp", rate = 1 / mean(x)), 1L),
      "gamma" = if (pos) {
        ft <- fd("gamma")
        if (!is.null(ft))
          add("gamma", list(family = "gamma",
                            shape = unname(ft$estimate["shape"]),
                            rate  = unname(ft$estimate["rate"])), 2L)
        else {
          m <- mean(x); v <- stats::var(x)
          if (v > 0)
            add("gamma", list(family = "gamma",
                              shape = m^2 / v, rate = m / v), 2L)
        }
      },
      "weibull" = if (pos) {
        ft <- fd("weibull")
        if (!is.null(ft))
          add("weibull", list(family = "weibull",
                              shape = unname(ft$estimate["shape"]),
                              scale = unname(ft$estimate["scale"])), 2L)
      },
      "logis" = {
        ft <- fd("logistic")
        if (!is.null(ft))
          add("logis", list(family = "logis",
                            location = unname(ft$estimate["location"]),
                            scale    = unname(ft$estimate["scale"])), 2L)
        else
          add("logis", list(family = "logis", location = mean(x),
                            scale = stats::sd(x) * sqrt(3) / pi), 2L)
      },
      "cauchy" = {
        ft <- fd("cauchy")
        if (!is.null(ft))
          add("cauchy", list(family = "cauchy",
                             location = unname(ft$estimate["location"]),
                             scale    = unname(ft$estimate["scale"])), 2L)
      },
      "pareto" = if (pos) {
        sc <- min(x)
        sh <- n / sum(log(x / sc))
        if (is.finite(sh) && sh > 0)
          add("pareto", list(family = "pareto",
                             shape = sh, scale = sc * (1 - 1e-9)), 2L)
      },
      "unif" = add("unif", list(family = "unif", min = a_j, max = b_j), 2L),
      "beta" = {
        rng <- b_j - a_j
        if (rng > 0) {
          u <- (x - a_j) / rng
          u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
          m <- mean(u); v <- stats::var(u)
          tt <- m * (1 - m) / v - 1
          if (is.finite(tt) && tt > 0)
            add("beta", list(family = "beta",
                             shape1 = m * tt, shape2 = (1 - m) * tt), 2L)
        }
      },
      # chisq / t / triang: no reliable generic estimator, user-specified only
      NULL
    )
  }
  out
}


# .fit_via_fitverse -----------------------------------------------------------
# Delegate the shared-parameterisation families to FitVerse and return its
# AIC-best model translated to the package's specification format.
.fit_via_fitverse <- function(x, families) {
  key_map <- c(norm = "normal",  lnorm = "lognormal", gamma = "gamma",
               exp  = "exponential", weibull = "weibull",
               cauchy = "cauchy", logis = "logistic", unif = "uniform")
  keys <- unname(key_map[intersect(families, names(key_map))])
  if (!length(keys)) return(NULL)

  ft <- tryCatch(
    FitVerse::fitverse(x, dists = keys, method = "MLE", criterion = "AIC",
                       plot = FALSE, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(ft) || is.null(ft$best_fit)) return(NULL)

  fam <- names(key_map)[match(ft$best_fit$dist_key, key_map)]
  if (is.na(fam)) return(NULL)
  pars <- as.list(ft$best_fit$params)

  spec <- c(list(family = fam), pars)
  ll   <- suppressWarnings(as.numeric(ft$best_fit$loglik))
  aic  <- suppressWarnings(as.numeric(ft$best_fit$aic))
  if (!is.finite(aic)) {
    # recompute on the same data with the package's own density
    pdf <- tryCatch(.build_pdf(spec, min(x), max(x)), error = function(e) NULL)
    if (is.null(pdf)) return(NULL)
    ll  <- sum(log(pmax(pdf(x), 1e-300)))
    aic <- 2 * length(pars) - 2 * ll
  }
  list(family = fam, spec = spec, ll = if (is.finite(ll)) ll else NA_real_,
       aic = aic, engine = "FitVerse")
}


# .fmt_params ------------------------------------------------------------------
.fmt_params <- function(spec) {
  pars <- spec[setdiff(names(spec), "family")]
  paste(sprintf("%s=%.4g", names(pars), unlist(pars)), collapse = ", ")
}
