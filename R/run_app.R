# =============================================================================
# run_app.R
#
# Launch the interactive Shiny application that ships with the package
# (installed under inst/shiny).  The user-interface packages (shiny, bslib, DT,
# jsonlite and the optional plotly / readxl / pagedown / FitVerse) are declared
# in 'Suggests', not 'Imports', so the core methodology package installs and
# passes R CMD check without them.  run_app() therefore uses them conditionally
# and reports clearly when a required one is missing.
# =============================================================================


#' Launch the multistrata Shiny Application
#'
#' Starts the interactive companion application bundled with the package.  The
#' app drives the whole design workflow without any coding: choose a data
#' source (a shipped example population, an uploaded CSV/Excel frame, a
#' simulated population, or a data-free \dQuote{densities only} specification),
#' pick the stratification and study variables, set the number of strata,
#' sample size or budget, priority weights and optional coefficient-of-variation
#' targets, and run the joint optimiser in empirical, parametric or
#' auto-parametric mode with the COBYLA or SLSQP solver.  Results appear as
#' stratum tables, allocation and boundary diagrams, interactive
#' visualisations, benchmark comparisons, precision (\code{find_n}) and
#' number-of-strata (\code{choose_H}) tools, an \eqn{\omega} trade-off explorer,
#' a downloadable report, and an auto-generated, reproducible R script.
#'
#' Every computation in the app is performed by this package's own exported
#' functions, so any design produced interactively can be reproduced exactly
#' from the console (the app's \dQuote{R Code} tab writes the script for you).
#'
#' @param launch.browser How to open the app.  A logical value is passed
#'   straight to \code{shiny::runApp()} (\code{TRUE} uses the default web
#'   browser, \code{FALSE} does not open one).  The string \code{"window"}
#'   (the default in interactive sessions) opens a borderless standalone
#'   Chromium/Chrome/Edge/Brave window when one is found and falls back to the
#'   default browser otherwise; \code{"browser"} forces the default browser.
#' @param port Optional integer port from 1 to 65535.  \code{NULL} (default)
#'   lets Shiny pick a free port.
#' @param ... Further arguments passed to \code{shiny::runApp()}.
#'
#' @return Invisibly \code{NULL}.  Called for its side effect: the function
#'   blocks the R session while the app is running and returns when it closes.
#'
#' @details
#' Required packages: \pkg{shiny}, \pkg{bslib}, \pkg{DT}, \pkg{jsonlite}.
#' Optional packages add features: \pkg{plotly} (interactive plots),
#' \pkg{readxl} (Excel upload), \pkg{pagedown} (direct PDF report export) and
#' \pkg{FitVerse} (an additional distribution-fitting engine).  Install the
#' required set once with
#' \code{install.packages(c("shiny", "bslib", "DT", "jsonlite"))}.
#'
#' @examples
#' \dontrun{
#'   run_app()                       # standalone window (or default browser)
#'   run_app(launch.browser = TRUE)  # force the default web browser
#'   run_app(port = 8123)
#' }
#'
#' @seealso \code{\link{stratify_multi}}, \code{\link{stratify_multi_density}},
#'   \code{\link{benchmark_designs}}
#'
#' @export
run_app <- function(launch.browser = c("window", "browser"),
                    port = NULL, ...) {

  required <- c("shiny", "bslib", "DT", "jsonlite")
  missing  <- required[!vapply(required, requireNamespace,
                               logical(1), quietly = TRUE)]
  if (length(missing))
    stop("The multistrata app requires the following package(s): ",
         paste(missing, collapse = ", "),
         ".\n  Install them with install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), ")).",
         call. = FALSE)

  app_dir <- system.file("shiny", package = "multistrata")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R")))
    stop("Could not locate the bundled Shiny app under inst/shiny. ",
         "Try reinstalling multistrata.", call. = FALSE)

  if (!is.null(port)) {
    port <- suppressWarnings(as.integer(port))
    if (is.na(port) || port < 1L || port > 65535L)
      stop("'port' must be an integer from 1 to 65535.", call. = FALSE)
  }

  # Resolve how the browser window is opened.
  if (missing(launch.browser)) {
    lb <- if (interactive()) .open_app_window
          else getOption("shiny.launch.browser", interactive())
  } else if (is.logical(launch.browser)) {
    lb <- isTRUE(launch.browser)
  } else {
    launch.browser <- match.arg(launch.browser)
    lb <- if (identical(launch.browser, "window")) .open_app_window else TRUE
  }

  shiny::runApp(app_dir, launch.browser = lb, port = port, ...)
  invisible(NULL)
}


# Open a borderless standalone window using a Chromium-based browser in
# "--app" mode; fall back to the default browser when none is found.
# Passed to shiny::runApp(launch.browser = ), so it receives the app URL.
.open_app_window <- function(url) {
  sysname  <- Sys.info()[["sysname"]]
  app_args <- c(paste0("--app=", url),
                "--window-size=1460,920",
                "--window-position=60,40")

  if (identical(sysname, "Darwin")) {
    browsers <- c(
      "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
      "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser",
      "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
      "/Applications/Chromium.app/Contents/MacOS/Chromium")
    for (b in browsers)
      if (file.exists(b)) {
        system2(b, args = app_args, wait = FALSE, stderr = FALSE)
        return(invisible(NULL))
      }
  } else if (identical(sysname, "Windows")) {
    candidates <- c(
      file.path(Sys.getenv("PROGRAMFILES"),
                "Google/Chrome/Application/chrome.exe"),
      file.path(Sys.getenv("PROGRAMFILES(X86)"),
                "Google/Chrome/Application/chrome.exe"),
      file.path(Sys.getenv("PROGRAMFILES(X86)"),
                "Microsoft/Edge/Application/msedge.exe"))
    for (b in candidates)
      if (file.exists(b)) {
        system2(b, args = app_args, wait = FALSE)
        return(invisible(NULL))
      }
  } else {
    for (b in c("google-chrome", "chromium-browser", "chromium", "brave-browser"))
      if (nzchar(Sys.which(b))) {
        system2(b, args = app_args, wait = FALSE)
        return(invisible(NULL))
      }
  }
  utils::browseURL(url)
}
