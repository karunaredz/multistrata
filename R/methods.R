# =============================================================================
# methods.R
#
# S3 methods for objects of class "multistrata".
#
# Fixes vs original:
#   print.multistrata(): boundary display used min/max of W_h as alpha/beta
#   (wrong values, wrong units).  Now uses alpha[j] and beta[j] from .meta
#   which are stored correctly in the fixed stratify_multi().
#   Variable labels now use actual column names via x_names.
# =============================================================================


#' Print method for multistrata
#' @param x     A \code{"multistrata"} object.
#' @param digits Decimal places for numeric output.  Default \code{4}.
#' @param ...   Unused.
#' @export
print.multistrata <- function(x, digits = 4, ...) {
  m <- x$.meta
  cat("\n")
  cat("================================================\n")
  cat("  Multivariate Optimum Stratification (multistrata)\n")
  cat("================================================\n")
  cat(sprintf("  Stratification variables : %d  (%s)\n",
              m$p, paste(m$x_names, collapse = ", ")))
  cat(sprintf("  Study variables          : %d  (%s)\n",
              m$q, paste(m$y_names, collapse = ", ")))
  cat(sprintf("  Strata per variable      : %s\n",
              paste(m$H_vec, collapse = " x ")))
  cat(sprintf("  Total strata             : %d\n",  m$H_total))
  cat(sprintf("  Population size (N)      : %d\n",  m$N))
  cat(sprintf("  Total sample (n)         : %d\n",  m$n_total))
  cat(sprintf("  Constraint type          : %s\n",  m$constraint))
  cat(sprintf("  Weight mode              : %s\n",
              if (!is.null(m$auto_fit)) "parametric (auto-selected)" else m$mode))
  if (!is.null(m$auto_fit) && !is.null(m$distributions))
    cat(sprintf("  Auto families            : %s\n",
        paste(sprintf("%s ~ %s", m$x_names,
                      vapply(m$distributions, function(d) d$family, "")),
              collapse = ", ")))
  cat(sprintf("  Boundary solver          : %s\n",
              toupper(if (is.null(m$method)) "cobyla" else m$method)))
  cat("------------------------------------------------\n")
  cat(sprintf("  Composite variance (Phi) : %.*f\n", digits, x$objective))
  cat(sprintf("  Optimality gap           : %.2f%%\n", 100 * x$optimality_gap))
  if (!is.null(x$kkt_residual))
    cat(sprintf("  KKT residual (scaled)    : %.2e  [%s]\n",
                x$kkt_residual,
                if (isTRUE(x$kkt_ok)) "< tol_kkt" else ">= tol_kkt"))
  if (!is.null(x$kkt_alloc) && !is.null(x$kkt_boundary))
    cat(sprintf("    allocation / boundary  : %.2e / %.2e\n",
                x$kkt_alloc, x$kkt_boundary))
  if (!is.null(x$cv_achieved)) {
    cat(sprintf("  CV achieved              : %s\n",
                paste(sprintf("%s = %.2f%%", m$y_names,
                              100 * x$cv_achieved), collapse = ", ")))
    if (!is.null(m$cv_targets))
      cat(sprintf("  CV targets               : %s\n",
                  paste(sprintf("%.2f%%", 100 * m$cv_targets), collapse = ", ")))
  }
  cat(sprintf("  Converged                : %s\n",
              if (x$converged) "YES" else "NO (max iterations reached)"))
  cat(sprintf("  Multi-starts used        : %d%s\n", x$n_starts_run,
              if (!is.null(x$n_starts_failed) && x$n_starts_failed > 0)
                sprintf(" (%d failed)", x$n_starts_failed) else ""))
  cat("------------------------------------------------\n")
  cat("  Optimal strata boundaries:\n")
  for (j in seq_len(m$p)) {
    bpts <- round(x$boundaries[[j]], digits)
    # FIX: use actual alpha[j] / beta[j] instead of min/max of W_h
    cat(sprintf("    %s: [%.4f | %s | %.4f]\n",
                m$x_names[j],
                m$alpha[j],
                paste(bpts, collapse = " | "),
                m$beta[j]))
  }
  cat("\n")
  invisible(x)
}


#' Summary method for multistrata
#' @param object A \code{"multistrata"} object.
#' @param ...    Unused.
#' @export
summary.multistrata <- function(object, ...) {
  m  <- object$.meta
  ss <- object$stratum_stats

  cat("\n")
  cat("== Stratum-Level Summary ========================\n")
  cat(sprintf("  omega (study-variable weights): %s\n",
              paste(round(m$omega, 3), collapse = ", ")))
  cat(sprintf("  Study variables: %s\n", paste(m$y_names, collapse = ", ")))
  cat("\n")

  disp <- data.frame(
    Stratum = ss$stratum,
    N_h     = ss$N_h,
    W_h     = round(ss$W_h, 4),
    n_h     = ss$n_h,
    A_h     = round(ss$A_h, 6),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(m$q)) {
    disp[[paste0("s2_", m$y_names[k])]] <- round(ss[[paste0("sigma2_k", k)]], 4)
  }

  print(disp, row.names = FALSE)
  cat(sprintf("\n  Sum n_h        = %d\n",   sum(ss$n_h)))
  cat(sprintf("  Composite Phi  = %.8f\n",   object$objective))
  cat(sprintf("  Optimality gap = %.2f%%\n", 100 * object$optimality_gap))
  cat("=================================================\n\n")
  invisible(object)
}


#' Basic plot method for multistrata (allocation or weights bar chart)
#'
#' For richer boundary visualisations use \code{\link{plot_strata}}.
#'
#' @param x     A \code{"multistrata"} object.
#' @param which One of \code{"allocation"} (default) or \code{"weights"}.
#' @param ...   Further graphical parameters passed to \code{barplot}.
#' @export
plot.multistrata <- function(x, which = "allocation", ...) {
  m  <- x$.meta
  ss <- x$stratum_stats

  # Validate before touching the graphics state. par() opens the default
  # device when none is active, so validating afterwards meant a mistyped
  # 'which' left an open device behind -- and in a non-interactive session
  # that is a stray Rplots.pdf in the working directory.
  if (!which %in% c("allocation", "weights"))
    stop("'which' must be 'allocation' or 'weights'.  ",
         "For boundary plots use plot_strata().")

  pal <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
           "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
  col_use <- pal[(seq_len(nrow(ss)) - 1L) %% length(pal) + 1L]

  op <- par(mar = c(7, 4.5, 3.5, 1.5))
  on.exit(par(op))

  if (which == "allocation") {
    barplot(height = ss$n_h, names.arg = ss$stratum, las = 2,
            col = col_use, border = "white",
            main = "Optimal Stratum Sample Allocation",
            ylab = expression(n[h] ~ "(sample size)"),
            cex.names = 0.62, ...)
    abline(h = mean(ss$n_h), col = "#E15759", lty = 2, lwd = 1.5)
    legend("topright", legend = sprintf("Mean = %.1f", mean(ss$n_h)),
           col = "#E15759", lty = 2, bty = "n", cex = 0.85)

  } else if (which == "weights") {
    barplot(height = ss$W_h, names.arg = ss$stratum, las = 2,
            col = col_use, border = "white",
            main = "Stratum Population Weights",
            ylab = expression(W[h] == N[h] / N),
            cex.names = 0.62, ...)
  } else {
    stop("'which' must be 'allocation' or 'weights'.  For boundary plots use plot_strata().")
  }
  invisible(x)
}
