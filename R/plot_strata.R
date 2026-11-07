# =============================================================================
# plot_strata.R
#
# Scientific visualisation of multivariate stratification boundaries produced
# by stratify_multi().
#
# Dispatch logic:
#   p == 1  -> .plot_strata_1d():
#              Kernel density + shaded strata + rug + per-stratum annotations.
#   p == 2  -> .plot_strata_2d():
#              TWO-PANEL figure.
#              Left, True 3-D scatter (base R persp + trans3d) with points
#                       coloured by stratum, semi-transparent boundary planes
#                       drawn at each b1/b2 cut-point, and stratum labels
#                       anchored at floor centroids.
#              Right, Clean 2-D top-down boundary grid with stratum name,
#                       n_h and W_h labelled inside each cell.
#   p >= 3  -> .plot_strata_nd():
#              p x p pairs matrix. Off-diagonal: scatter + boundary lines.
#              Diagonal: kernel density histogram.
#
# Dependencies: base R only (stats, graphics, grDevices).
# =============================================================================


#' Visualise Multivariate Stratification Boundaries
#'
#' Produces a scientific, publication-quality diagram of the optimal strata
#' boundaries returned by \code{\link{stratify_multi}}.
#'
#' \strong{p = 1}: Kernel density with shaded stratum regions, boundary
#' markers, and within-stratum statistics.\cr
#' \strong{p = 2}: Two-panel figure.
#' \emph{Left} - A true 3-D scatter plot (base R perspective projection) of
#' population units coloured by stratum assignment.  Semi-transparent planar
#' walls are drawn at each boundary cut-point, slicing the data cloud into
#' labelled stratum regions.  The vertical (Z) axis is the first study
#' variable.
#' \emph{Right} - A clean 2-D top-down boundary grid for reference, with each
#' cell labelled by stratum name, sample allocation \eqn{n_h}, and population
#' weight \eqn{W_h}.\cr
#' \strong{p >= 3}: \eqn{p \times p} pairs matrix.  Off-diagonal panels show
#' 2-D projections of boundaries onto the relevant variable pair.  Diagonal
#' panels show kernel density histograms.
#'
#' @param result     A \code{"multistrata"} object from
#'   \code{\link{stratify_multi}}.
#' @param data       The original stratification variable matrix/data.frame
#'   (same object passed to \code{stratify_multi}).
#' @param study_vars The original study-variable matrix (same object passed
#'   to \code{stratify_multi}).  Required when \code{p = 2} (used as the
#'   Z-axis of the 3-D scatter).  May be \code{NULL} only when \code{p != 2}.
#' @param jitter_n   Maximum scatter points displayed.  A random sub-sample
#'   is drawn when \code{nrow(data) > jitter_n}.  Default \code{3000}.
#' @param alpha_pt   Point transparency for scatter plots.  Default \code{0.40}.
#' @param pt_cex     Point size.  Default \code{0.7}.
#' @param palette    Character vector of colours, one per stratum (recycled).
#'   Defaults to a high-contrast qualitative palette.
#' @param theta_3d   Azimuthal viewing angle for the 3-D scatter (\code{p=2}).
#'   Default \code{220}.
#' @param phi_3d     Colatitude viewing angle.  Default \code{22}.
#' @param show_allocation Logical.  Print \eqn{n_h}/\eqn{W_h} in the 2-D
#'   reference panel (\code{p=2}) and in the 1-D annotations.  Default
#'   \code{TRUE}.
#' @param ...        Additional graphical parameters passed to the primary
#'   \code{persp} or \code{plot} call.
#'
#' @return Invisibly returns \code{result}.
#'
#' @examples
#' set.seed(42)
#' pop <- sim_population(N = 3000, p = 2, q = 2, scenario = "independent")
#' res <- stratify_multi(
#'   data       = pop[, c("X1","X2")],
#'   study_vars = pop[, c("Y1","Y2")],
#'   H_vec      = c(3, 3), n_total = 300, n_starts = 5
#' )
#' plot_strata(res, data = pop[, c("X1","X2")],
#'             study_vars = pop[, c("Y1","Y2")])
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{plot.multistrata}}
#' @export
plot_strata <- function(result,
                        data,
                        study_vars      = NULL,
                        jitter_n        = 3000L,
                        alpha_pt        = 0.40,
                        pt_cex          = 0.70,
                        palette         = NULL,
                        theta_3d        = 220,
                        phi_3d          = 22,
                        show_allocation = TRUE,
                        ...) {

  if (!inherits(result, "multistrata"))
    stop("'result' must be a 'multistrata' object from stratify_multi().")

  data <- as.matrix(data)
  if (!is.null(study_vars)) study_vars <- as.matrix(study_vars)

  m <- result$.meta
  p <- m$p

  # (study_vars is accepted but not required for p == 2 with a density surface)

  H_total  <- m$H_total
  ss       <- result$stratum_stats
  strat_id <- .assign_strata(data, result$boundaries, m$H_vec)

  # Colour palette - one per stratum
  def_pal <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
    "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
  )
  pal        <- if (is.null(palette)) def_pal else palette
  col_strata <- pal[(seq_len(H_total) - 1L) %% length(pal) + 1L]

  if (p == 1L) {
    .plot_strata_1d(result, data, strat_id, ss, m, col_strata, show_allocation)
  } else if (p == 2L) {
    .plot_strata_2d(result, data, study_vars, strat_id, ss, m,
                    col_strata, jitter_n, alpha_pt, pt_cex,
                    theta_3d, phi_3d, show_allocation, ...)
  } else {
    .plot_strata_nd(result, data, strat_id, ss, m,
                    col_strata, jitter_n, alpha_pt, pt_cex, ...)
  }
  invisible(result)
}


# =============================================================================
# 1-D renderer
# =============================================================================
.plot_strata_1d <- function(result, data, strat_id, ss, m, col_strata,
                             show_allocation) {
  x      <- data[, 1L]
  B      <- result$boundaries[[1L]]
  H      <- m$H_vec[1L]
  breaks <- c(m$alpha[1L], B, m$beta[1L])

  op <- par(mar = c(5, 4.5, 4.5, 2))
  on.exit(par(op))

  dens <- density(x, n = 1024L)
  ylim <- c(0, max(dens$y) * 1.18)

  plot(dens, main = "", xlab = m$x_names[1L], ylab = "Density",
       col = "grey30", lwd = 2, bty = "l", ylim = ylim)

  # Shaded stratum regions
  for (h in seq_len(H)) {
    lo <- breaks[h];  hi <- breaks[h + 1L]
    xi <- dens$x[dens$x >= lo & dens$x <= hi]
    yi <- dens$y[dens$x >= lo & dens$x <= hi]
    if (length(xi) > 1L)
      polygon(c(lo, xi, hi), c(0, yi, 0),
              col = adjustcolor(col_strata[h], alpha.f = 0.32), border = NA)
  }
  lines(dens, col = "grey20", lwd = 2)

  # Boundary lines + top-axis labels
  abline(v = B, col = "grey10", lwd = 2.2)
  for (b in B)
    mtext(sprintf("%.3f", b), side = 3, at = b, cex = 0.72,
          col = "grey10", line = 0.3)

  # Coloured rug
  for (h in seq_len(H))
    rug(x[strat_id == h],
        col = adjustcolor(col_strata[h], alpha.f = 0.55),
        ticksize = 0.014, lwd = 0.5)

  # Per-stratum annotations
  for (h in seq_len(H)) {
    mid   <- (breaks[h] + breaks[h + 1L]) / 2
    y_ann <- max(dens$y) * (0.95 - 0.14 * ((h - 1L) %% 2L))
    lab   <- if (show_allocation)
      sprintf("S%d\nn_h=%d\nsigma2=%.2f", h, ss$n_h[h], ss$sigma2_k1[h])
    else
      sprintf("S%d", h)
    text(mid, y_ann, labels = lab, cex = 0.72, col = col_strata[h], font = 2)
  }

  legend("topright",
         legend = paste0("Stratum ", seq_len(H)),
         fill   = adjustcolor(col_strata[seq_len(H)], alpha.f = 0.55),
         border = NA, bty = "n", cex = 0.82)

  title(main = sprintf(
    "Optimal Stratification  |  %s  |  H = %d  |  Phi* = %.6f",
    m$x_names[1L], H, result$objective), cex.main = 0.92)
}


# =============================================================================
# 2-D renderer
# Left panel : 3-D joint density surface (MASS::kde2d) coloured by stratum,
#              with semi-transparent boundary planes at each cut value and
#              stratum labels at floor centroids.
# Right panel: 2-D top-down boundary grid reference with stratum name, n_h, W_h.
# =============================================================================
.plot_strata_2d <- function(result, data, study_vars, strat_id, ss, m,
                             col_strata, jitter_n, alpha_pt, pt_cex,
                             theta_3d, phi_3d, show_allocation, ...) {

  B1 <- result$boundaries[[1L]];  B2 <- result$boundaries[[2L]]
  H1 <- m$H_vec[1L];  H2 <- m$H_vec[2L]
  N  <- nrow(data)

  # Axis ranges with 5% padding, and full boundary vectors
  x1r <- range(data[, 1L]);  x2r <- range(data[, 2L])
  x1p <- diff(x1r) * 0.05;   x2p <- diff(x2r) * 0.05
  x1_lim <- x1r + c(-x1p,  x1p)
  x2_lim <- x2r + c(-x2p,  x2p)
  bx <- c(m$alpha[1L], B1, m$beta[1L])
  by <- c(m$alpha[2L], B2, m$beta[2L])

  # Left panel wider (3:2 ratio)
  op <- par(oma = c(0, 0, 2.8, 0))
  on.exit(par(op))
  layout(matrix(c(1L, 2L), 1L, 2L), widths = c(3, 2))

  # ============================================================
  # Panel 1: 3-D joint density surface
  # ============================================================
  par(mar = c(2, 1, 3.5, 1))

  # 2-D kernel density estimate on a fine grid
  grid_n  <- 60L
  dens2d  <- MASS::kde2d(data[, 1L], data[, 2L], n = grid_n,
                          lims = c(x1_lim, x2_lim))
  gx   <- dens2d$x   # length grid_n
  gy   <- dens2d$y   # length grid_n
  gz   <- dens2d$z   # grid_n x grid_n
  z_max <- max(gz) * 1.05   # ceiling for boundary walls

  # Colour each surface facet by the stratum its centre falls in.
  # Stratum index: 1 + idx1 + idx2 * H1  (matches .assign_strata() mult formula)
  facet_col <- matrix(NA_character_, grid_n - 1L, grid_n - 1L)
  for (i in seq_len(grid_n - 1L)) {
    cx1  <- (gx[i] + gx[i + 1L]) / 2
    idx1 <- findInterval(cx1, B1)          # 0 .. H1-1
    for (j in seq_len(grid_n - 1L)) {
      cx2  <- (gy[j] + gy[j + 1L]) / 2
      idx2 <- findInterval(cx2, B2)        # 0 .. H2-1
      h    <- 1L + idx1 + idx2 * H1
      facet_col[i, j] <- adjustcolor(col_strata[h], alpha.f = 0.80)
    }
  }

  # Draw the 3-D density surface
  pmat <- persp(
    x        = gx,
    y        = gy,
    z        = gz,
    col      = facet_col,
    border   = NA,
    shade    = 0.20,
    theta    = theta_3d,
    phi      = phi_3d,
    expand   = 0.55,
    xlab     = m$x_names[1L],
    ylab     = m$x_names[2L],
    zlab     = "Density",
    ticktype = "detailed",
    cex.axis = 0.58,
    cex.lab  = 0.80,
    main     = "",
    zlim     = c(0, z_max),
    ...
  )

  # Helper: project and fill a 3-D polygon
  .poly3d <- function(px1, px2, pz, col, border = NA, lwd = 1.5, lty = 1) {
    pp <- trans3d(px1, px2, pz, pmat)
    polygon(pp$x, pp$y, col = col, border = border, lwd = lwd, lty = lty)
  }

  # -- Boundary walls perpendicular to X1 axis (solid planes) -------------
  wx2 <- c(x2_lim[1L], x2_lim[2L], x2_lim[2L], x2_lim[1L])
  wz  <- c(0,          0,           z_max,       z_max)
  for (b in B1) {
    .poly3d(rep(b, 4L), wx2, wz,
            col    = adjustcolor("grey30", alpha.f = 0.14),
            border = adjustcolor("grey10", alpha.f = 0.65),
            lwd    = 2.0)
    lp  <- trans3d(b, x2_lim[1L], z_max, pmat)
    lab <- sprintf("%s = %.3f", m$x_names[1L], b)
    lw  <- strwidth(lab, cex = 0.62, font = 2)
    lh  <- strheight(lab, cex = 0.62, font = 2)
    rect(lp$x - 0.55 * lw, lp$y + 0.15 * lh, lp$x + 0.55 * lw,
         lp$y + 1.55 * lh, col = adjustcolor("white", alpha.f = 0.75),
         border = NA, xpd = NA)
    text(lp$x, lp$y, labels = lab,
         cex = 0.62, col = "grey10", font = 2, adj = c(0.5, -0.35), xpd = NA)
  }

  # -- Boundary walls perpendicular to X2 axis (dashed planes) ------------
  wx1 <- c(x1_lim[1L], x1_lim[2L], x1_lim[2L], x1_lim[1L])
  for (b in B2) {
    .poly3d(wx1, rep(b, 4L), wz,
            col    = adjustcolor("grey30", alpha.f = 0.11),
            border = adjustcolor("grey20", alpha.f = 0.55),
            lwd    = 2.0, lty = 2L)
    lp  <- trans3d(x1_lim[2L], b, z_max, pmat)
    lab <- sprintf("%s = %.3f", m$x_names[2L], b)
    lw  <- strwidth(lab, cex = 0.62, font = 2)
    lh  <- strheight(lab, cex = 0.62, font = 2)
    rect(lp$x - 0.55 * lw, lp$y + 0.15 * lh, lp$x + 0.55 * lw,
         lp$y + 1.55 * lh, col = adjustcolor("white", alpha.f = 0.75),
         border = NA, xpd = NA)
    text(lp$x, lp$y, labels = lab,
         cex = 0.62, col = "grey20", font = 2, adj = c(0.5, -0.35), xpd = NA)
  }

  # -- Stratum labels at floor (z = 0) centroid of each region ------------
  # Drawn on a translucent white halo so they remain readable over the
  # shaded floor tiles and boundary walls.
  h_cnt <- 0L
  for (i in seq_len(H1)) {
    for (j in seq_len(H2)) {
      h_cnt <- h_cnt + 1L
      cx <- (bx[i] + bx[i + 1L]) / 2
      cy <- (by[j] + by[j + 1L]) / 2
      lp <- trans3d(cx, cy, 0, pmat)
      lab <- sprintf("%s\nn=%d", ss$stratum[h_cnt], ss$n_h[h_cnt])
      lw <- strwidth(lab,  cex = 0.58, font = 2)
      lh <- strheight(lab, cex = 0.58, font = 2)
      rect(lp$x - 0.55 * lw, lp$y - 0.60 * lh,
           lp$x + 0.55 * lw, lp$y + 0.60 * lh,
           col = adjustcolor("white", alpha.f = 0.60), border = NA, xpd = NA)
      text(lp$x, lp$y, labels = lab,
           cex = 0.58, font = 2, col = col_strata[h_cnt],
           adj = c(0.5, 0.5), xpd = NA)
    }
  }

  # -- Legend ---------------------------------------------------------------
  legend("topleft",
         legend = ss$stratum,
         fill   = adjustcolor(col_strata[seq_len(m$H_total)], alpha.f = 0.80),
         border = NA, bty = "n", cex = 0.62,
         ncol   = max(1L, ceiling(m$H_total / 5L)))

  title(main = sprintf("Joint density: f(%s, %s)  |  strata coloured",
                       m$x_names[1L], m$x_names[2L]),
        cex.main = 0.88, col.main = "grey20")

  # ============================================================
  # Panel 2: 2-D top-down boundary grid (reference)
  # ============================================================
  par(mar = c(3, 3, 3.5, 2))

  # Sub-sample for scatter overlay
  idx_plot <- if (N > jitter_n) sample.int(N, jitter_n) else seq_len(N)

# FIX: axis limits must span the full stratum grid (alpha..beta), not just
  # the sub-sampled scatter, otherwise the centres of large outer cells fall
  # on or beyond the plot edge and their labels are clipped.
  plot(data[idx_plot, 1L], data[idx_plot, 2L],
       col  = adjustcolor(col_strata[strat_id[idx_plot]], alpha.f = alpha_pt * 0.6),
       pch  = 16, cex = pt_cex * 0.75,
       xlim = range(bx), ylim = range(by),
       xlab = m$x_names[1L], ylab = m$x_names[2L],
       main = "Boundary grid (top view)",
       bty  = "l", cex.main = 0.88)

  # Shaded cells
  h_cnt <- 0L
  for (i in seq_len(H1)) {
    for (j in seq_len(H2)) {
      h_cnt <- h_cnt + 1L
      rect(bx[i], by[j], bx[i + 1L], by[j + 1L],
           col = adjustcolor(col_strata[h_cnt], alpha.f = 0.18), border = NA)
    }
  }

  # Boundary lines
  abline(v = B1, col = "grey15", lwd = 2)
  abline(h = B2, col = "grey15", lwd = 2)

  # Boundary value labels on outer axes
  axis(3, at = B1, labels = sprintf("%.3f", B1),
       cex.axis = 0.68, col.axis = "grey10", tick = FALSE, line = -0.5)
  axis(4, at = B2, labels = sprintf("%.3f", B2),
       cex.axis = 0.68, col.axis = "grey10", tick = FALSE, line = -0.5)

  # Per-cell labels: stratum name + n_h + W_h, centred in each cell,
  # clamped inside the plot region, on a translucent white halo box so
  # they stay readable over shading and scatter points.
  usr <- par("usr")
  h_cnt <- 0L
  for (i in seq_len(H1)) {
    for (j in seq_len(H2)) {
      h_cnt <- h_cnt + 1L
      cx <- (bx[i] + bx[i + 1L]) / 2
      cy <- (by[j] + by[j + 1L]) / 2
      lab <- if (show_allocation)
        sprintf("%s\nn=%d\nW=%.3f", ss$stratum[h_cnt], ss$n_h[h_cnt], ss$W_h[h_cnt])
      else
        ss$stratum[h_cnt]
      lw <- strwidth(lab,  cex = 0.68, font = 2)
      lh <- strheight(lab, cex = 0.68, font = 2)
      # clamp so the whole label box stays inside the plotting region
      cx <- min(max(cx, usr[1L] + 0.55 * lw), usr[2L] - 0.55 * lw)
      cy <- min(max(cy, usr[3L] + 0.60 * lh), usr[4L] - 0.60 * lh)
      rect(cx - 0.55 * lw, cy - 0.60 * lh, cx + 0.55 * lw, cy + 0.60 * lh,
           col = adjustcolor("white", alpha.f = 0.70), border = NA)
      text(cx, cy, labels = lab, cex = 0.68, font = 2,
           col = col_strata[h_cnt], xpd = NA)
    }
  }

  # Overall figure title
  title(outer = TRUE,
        main  = sprintf(
          "multistrata  |  H = %s  |  n = %d  |  Phi* = %.6f  |  gap = %.2f%%",
          paste(m$H_vec, collapse = "x"), m$n_total,
          result$objective, 100 * result$optimality_gap),
        cex.main = 0.90, col.main = "grey20")
}


# =============================================================================
# n-D renderer  (p >= 3) -- p x p pairs matrix
# Off-diagonal: scatter coloured by stratum + boundary lines for each axis var
# Diagonal:     kernel density + histogram overlay
# =============================================================================
.plot_strata_nd <- function(result, data, strat_id, ss, m,
                             col_strata, jitter_n, alpha_pt, pt_cex, ...) {
  p <- m$p
  N <- nrow(data)

  idx_plot <- if (N > jitter_n) sample.int(N, jitter_n) else seq_len(N)
  d_plot   <- data[idx_plot, , drop = FALSE]
  sid_plot <- strat_id[idx_plot]
  pt_col   <- adjustcolor(col_strata[sid_plot], alpha.f = alpha_pt)

  op <- par(mfrow = c(p, p), mar = c(2, 2, 1, 0.5), oma = c(0, 0, 3, 0))
  on.exit(par(op))

  for (row in seq_len(p)) {
    for (col in seq_len(p)) {

      if (row == col) {
        # Diagonal: kernel density histogram
        x_d  <- data[, row]
        dens <- density(x_d, n = 512L)
        plot(dens, main = "", xlab = "", ylab = "", bty = "l",
             col  = col_strata[1L], lwd = 1.8,
             yaxt = "n", cex.axis = 0.65)
        h_d <- hist(x_d, plot = FALSE, breaks = "FD")
        h_d$density <- h_d$density / max(h_d$density) * max(dens$y)
        rect(h_d$breaks[-length(h_d$breaks)], 0,
             h_d$breaks[-1L], h_d$density,
             col = adjustcolor(col_strata[1L], 0.25), border = NA)
        lines(dens, col = col_strata[1L], lwd = 1.8)
        mtext(m$x_names[row], side = 3, cex = 0.75, font = 2, line = 0)

      } else {
        # Off-diagonal: scatter + boundary lines for each axis variable
        plot(d_plot[, col], d_plot[, row],
             col  = pt_col, pch = 16, cex = pt_cex,
             xlab = "", ylab = "", bty = "l",
             cex.axis = 0.60, ...)

        # Faint cell shading for this 2-D projection
        Hx <- m$H_vec[col];  Hy <- m$H_vec[row]
        bx <- c(m$alpha[col], result$boundaries[[col]], m$beta[col])
        by <- c(m$alpha[row], result$boundaries[[row]], m$beta[row])
        h_cnt <- 0L
        for (hx in seq_len(Hx)) {
          for (hy in seq_len(Hy)) {
            h_cnt <- h_cnt + 1L
            rect(bx[hx], by[hy], bx[hx + 1L], by[hy + 1L],
                 col    = adjustcolor(
                   col_strata[(h_cnt - 1L) %% m$H_total + 1L], alpha.f = 0.07),
                 border = NA)
          }
        }

        # Boundary lines
        if (length(result$boundaries[[col]]) > 0L)
          abline(v = result$boundaries[[col]], col = "grey10", lwd = 1.6)
        if (length(result$boundaries[[row]]) > 0L)
          abline(h = result$boundaries[[row]], col = "grey10", lwd = 1.6)
      }
    }
  }

  # Overall title
  title(outer = TRUE,
        main  = sprintf(
          "multistrata  |  H = %s  |  n = %d  |  Phi* = %.6f  |  gap = %.2f%%",
          paste(m$H_vec, collapse = "x"), m$n_total,
          result$objective, 100 * result$optimality_gap),
        cex.main = 0.92, col.main = "grey20")

  # Compact stratum legend
  if (m$H_total <= 12L) {
    op2 <- par(fig = c(0.85, 1, 0.85, 1), new = TRUE, mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
           legend = ss$stratum,
           col    = col_strata[seq_len(m$H_total)],
           pch    = 16, pt.cex = 1.1, bty = "n", cex = 0.65,
           ncol   = max(1L, ceiling(m$H_total / 6L)),
           title  = "Strata")
    par(op2)
  }
}
