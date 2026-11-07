# =============================================================================
# plot_strata_rgl.R
#
# Interactive 3-D density surface for a bivariate (p = 2) stratification
# result.  Requires the 'rgl' package (listed in Suggests).
#
# The function opens a live rgl window.  Use the mouse to:
#   Left-drag   : rotate
#   Right-drag  : zoom
#   Middle-drag : pan
#   Scroll      : zoom
#
# Rendering:
#   - Joint density surface f(X1, X2) from MASS::kde2d, facets coloured
#     by stratum assignment.
#   - Semi-transparent planar walls at each boundary cut-point.
#     Walls for X1 boundaries are solid; walls for X2 boundaries are
#     drawn with a slightly different colour so they are distinguishable.
#   - Stratum name + n_h label at the floor centroid of each region.
#   - Axes, title, and a 2-D legend inset.
# =============================================================================


#' Interactive 3-D Density Surface for Multivariate Stratification
#'
#' Opens a rotatable \pkg{rgl} window showing the joint density surface
#' \eqn{f(X_1, X_2)} coloured by stratum, with semi-transparent boundary
#' planes at each optimal cut-point and stratum labels at floor centroids.
#'
#' Use the mouse to rotate (left-drag), zoom (scroll or right-drag), and
#' pan (middle-drag).
#'
#' @param result         A \code{"multistrata"} object (\code{p = 2} only).
#' @param data           Stratification variable matrix/data.frame (same
#'   object passed to \code{\link{stratify_multi}}).
#' @param palette        Optional colour vector, one colour per stratum
#'   (recycled).  Defaults to the same high-contrast palette used by
#'   \code{\link{plot_strata}}.
#' @param grid_n         KDE grid resolution.  Default \code{60}.
#' @param alpha_surface  Opacity of the density surface (0 = transparent,
#'   1 = opaque).  Default \code{0.82}.
#' @param alpha_plane    Opacity of the boundary plane walls.  Default
#'   \code{0.22}.
#' @param show_allocation Logical.  Include \eqn{n_h} in the stratum floor
#'   labels.  Default \code{TRUE}.
#' @param windowTitle    Title shown in the rgl window's title bar.
#' @param ...            Additional arguments forwarded to
#'   \code{rgl::persp3d()}.
#'
#' @return Invisibly returns \code{result}.  The rgl window remains open
#'   for interaction.
#'
#' @examples
#' \dontrun{
#' set.seed(2025)
#' pop <- sim_population(N = 5000, p = 2, q = 2,
#'                        scenario = "correlated", rho = 0.6)
#' res <- stratify_multi(
#'   data = pop[, c("X1","X2")], study_vars = pop[, c("Y1","Y2")],
#'   H_vec = c(3, 3), n_total = 500, n_starts = 10)
#' plot_strata_rgl(res, data = pop[, c("X1","X2")])
#' }
#'
#' @seealso \code{\link{plot_strata}}, \code{\link{stratify_multi}}
#' @export
plot_strata_rgl <- function(result,
                             data,
                             palette          = NULL,
                             grid_n           = 60L,
                             alpha_surface    = 0.82,
                             alpha_plane      = 0.22,
                             show_allocation  = TRUE,
                             windowTitle      = "multistrata - interactive",
                             ...) {

  # ── dependency check ────────────────────────────────────────────────────────
  # Distinguish "rgl is absent" from "rgl is installed but will not load".
  # The second case is common on recent macOS, which no longer ships the
  # OpenGL libraries rgl links against, so loading fails with a dlopen error
  # on libGLU. Telling such a user to install rgl is useless advice -- they
  # already have it -- so point them at the static 3-D view instead, which
  # needs no OpenGL. Loading rgl emits its own warnings on the way to
  # failing; they are captured here so the message the user sees is ours.
  if (!.rgl_available()) {
    rgl_installed <- .rgl_installed()
    stop(if (rgl_installed)
      paste0("The 'rgl' package is installed but cannot be loaded, so the ",
             "interactive 3-D view is unavailable.\n",
             "  This is almost always a missing OpenGL library (rgl needs ",
             "libGLU); recent macOS releases no longer provide one.\n",
             "  Use plot_strata(result, data, study_vars) instead: for two ",
             "stratification variables it draws a static 3-D perspective ",
             "plot with the boundary walls, and needs no OpenGL.")
      else
      paste0("Package 'rgl' is required for the interactive 3-D view.\n",
             "  Install it with:  install.packages(\"rgl\")\n",
             "  Or use plot_strata(result, data, study_vars) for the static ",
             "3-D perspective plot, which needs no extra packages."),
      call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Package 'MASS' is required for kde2d.\n",
         "  Install it with:  install.packages('MASS')")

  if (!inherits(result, "multistrata"))
    stop("'result' must be a 'multistrata' object from stratify_multi().")

  data <- as.matrix(data)
  m    <- result$.meta

  if (m$p != 2L)
    stop("plot_strata_rgl() currently supports p = 2 stratification ",
         "variables only.")

  B1 <- result$boundaries[[1L]];  B2 <- result$boundaries[[2L]]
  H1 <- m$H_vec[1L];  H2 <- m$H_vec[2L]
  ss <- result$stratum_stats

  # ── colour palette ───────────────────────────────────────────────────────────
  def_pal <- c(
    "#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
    "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD",
    "#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF"
  )
  pal        <- if (is.null(palette)) def_pal else palette
  col_strata <- pal[(seq_len(m$H_total) - 1L) %% length(pal) + 1L]

  # ── axis limits and boundary vectors ────────────────────────────────────────
  x1r <- range(data[, 1L]);  x2r <- range(data[, 2L])
  x1p <- diff(x1r) * 0.05;  x2p <- diff(x2r) * 0.05
  x1_lim <- x1r + c(-x1p,  x1p)
  x2_lim <- x2r + c(-x2p,  x2p)
  bx <- c(m$alpha[1L], B1, m$beta[1L])
  by <- c(m$alpha[2L], B2, m$beta[2L])

  # ── 2-D KDE ─────────────────────────────────────────────────────────────────
  dens2d <- MASS::kde2d(data[, 1L], data[, 2L],
                         n    = grid_n,
                         lims = c(x1_lim, x2_lim))
  gx    <- dens2d$x
  gy    <- dens2d$y
  gz    <- dens2d$z
  z_max <- max(gz) * 1.10   # ceiling height for boundary walls

  # ── facet colours: one colour per (grid_n-1) x (grid_n-1) surface patch ─────
  # Stratum index: 1 + idx1 + idx2 * H1  (matches .assign_strata mult formula)
  facet_col <- matrix(NA_character_, grid_n - 1L, grid_n - 1L)
  for (i in seq_len(grid_n - 1L)) {
    cx1  <- (gx[i] + gx[i + 1L]) / 2
    idx1 <- findInterval(cx1, B1)        # 0 .. H1-1
    for (j in seq_len(grid_n - 1L)) {
      cx2  <- (gy[j] + gy[j + 1L]) / 2
      idx2 <- findInterval(cx2, B2)      # 0 .. H2-1
      h    <- 1L + idx1 + idx2 * H1
      facet_col[i, j] <- col_strata[h]
    }
  }

  # ── open rgl window ──────────────────────────────────────────────────────────
  rgl::open3d(windowTitle = windowTitle)
  rgl::bg3d(color = "white")

  # ── density surface ──────────────────────────────────────────────────────────
  rgl::persp3d(
    x      = gx,
    y      = gy,
    z      = gz,
    col    = facet_col,
    alpha  = alpha_surface,
    xlab   = m$x_names[1L],
    ylab   = m$x_names[2L],
    zlab   = "Density",
    smooth = TRUE,
    ...
  )

  # ── boundary planes at X1 cut-points (solid colour) ─────────────────────────
  for (b in B1) {
    # Four corners of the vertical wall: (b, y_min/max, 0/z_max)
    rgl::quads3d(
      x   = rep(b, 4L),
      y   = c(x2_lim[1L], x2_lim[2L], x2_lim[2L], x2_lim[1L]),
      z   = c(0,           0,           z_max,       z_max),
      col = "#2c2c2c",
      alpha = alpha_plane,
      lit = FALSE
    )
    # Label above the wall
    rgl::text3d(
      x     = b,
      y     = x2_lim[1L],
      z     = z_max * 1.06,
      texts = sprintf("%s = %.3f", m$x_names[1L], b),
      col   = "grey10",
      cex   = 0.88,
      font  = 2
    )
  }

  # ── boundary planes at X2 cut-points (lighter colour, dashed feel) ───────────
  for (b in B2) {
    rgl::quads3d(
      x   = c(x1_lim[1L], x1_lim[2L], x1_lim[2L], x1_lim[1L]),
      y   = rep(b, 4L),
      z   = c(0,           0,           z_max,       z_max),
      col = "#555555",
      alpha = alpha_plane * 0.80,
      lit = FALSE
    )
    rgl::text3d(
      x     = x1_lim[2L],
      y     = b,
      z     = z_max * 1.06,
      texts = sprintf("%s = %.3f", m$x_names[2L], b),
      col   = "grey25",
      cex   = 0.88,
      font  = 2
    )
  }

  # ── stratum floor labels (centroid of each region at z = 0) ─────────────────
  h_cnt <- 0L
  for (i in seq_len(H1)) {
    for (j in seq_len(H2)) {
      h_cnt <- h_cnt + 1L
      cx <- (bx[i] + bx[i + 1L]) / 2
      cy <- (by[j] + by[j + 1L]) / 2
      lab <- if (show_allocation)
        paste0(ss$stratum[h_cnt], "\nn=", ss$n_h[h_cnt])
      else
        ss$stratum[h_cnt]
      rgl::text3d(
        x = cx, y = cy, z = 0,
        texts = lab,
        col   = col_strata[h_cnt],
        cex   = 0.82,
        font  = 2
      )
    }
  }

  # ── axes ─────────────────────────────────────────────────────────────────────
  rgl::axes3d(edges = c("x--", "y--", "z--"), col = "grey40", nticks = 5)
  rgl::box3d(col = "grey80")

  # ── title ────────────────────────────────────────────────────────────────────
  rgl::title3d(
    main = sprintf(
      "H = %s  |  n = %d  |  \u03a6* = %.5f  |  gap = %.2f%%",
      paste(m$H_vec, collapse = " x "),
      m$n_total,
      result$objective,
      100 * result$optimality_gap),
    xlab = m$x_names[1L],
    ylab = m$x_names[2L],
    zlab = "Density",
    col  = "grey15",
    cex  = 0.82
  )

  # ── 2-D legend inset ─────────────────────────────────────────────────────────
  rgl::legend3d(
    "topright",
    legend  = ss$stratum,
    fill    = col_strata[seq_len(m$H_total)],
    border  = NA,
    bty     = "n",
    cex     = 0.72,
    title   = "Strata",
    ncol    = max(1L, ceiling(m$H_total / 5L))
  )

  invisible(result)
}


# rgl availability --------------------------------------------------------
# Split into two tiny helpers so the failure branches of plot_strata_rgl()
# are testable. The interesting case -- rgl present on disk but unable to
# load -- cannot be reproduced on a machine that has working OpenGL, so the
# tests mock these rather than the platform.

# TRUE when rgl can actually be loaded and used. requireNamespace() runs
# rgl's .onLoad, which on a machine without OpenGL emits warnings before
# returning FALSE; they are suppressed so the caller sees our message and
# not rgl's dlopen trace.
.rgl_available <- function()
  isTRUE(suppressWarnings(suppressMessages(
    requireNamespace("rgl", quietly = TRUE))))

# TRUE when rgl is installed, whether or not it can be loaded. Deliberately
# does not attempt to load it.
.rgl_installed <- function() nzchar(system.file(package = "rgl"))
