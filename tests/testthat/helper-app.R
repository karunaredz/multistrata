# Helpers for exercising the bundled Shiny app (inst/shiny/).
# The UI packages live in Suggests, so every app test skips gracefully when
# they (or the installed app files) are unavailable.

skip_if_no_app <- function() {
  testthat::skip_if_not_installed("shiny")
  testthat::skip_if_not_installed("bslib")
  testthat::skip_if_not_installed("DT")
  if (!nzchar(system.file("shiny", "app.R", package = "multistrata")))
    testthat::skip("bundled Shiny app is not installed")
}

# Source the bundled app into a fresh environment (defines ui, server and the
# app_helpers.R functions via source(local = TRUE)); returns that environment.
#
# `without` names optional packages that the app should believe are absent. It
# works by shadowing requireNamespace() inside the app's environment, which is
# how the app probes for its Suggests, so the has_* flags come out FALSE
# without the package actually having to be uninstalled. This lets the suite
# exercise the degraded paths on a fully provisioned machine.
load_app_env <- function(without = character()) {
  f   <- system.file("shiny", "app.R", package = "multistrata")
  env <- new.env(parent = globalenv())
  if (length(without))
    env$requireNamespace <- function(package, ..., quietly = FALSE)
      if (package %in% without) FALSE
      else base::requireNamespace(package, ..., quietly = quietly)
  old <- setwd(dirname(f))
  on.exit(setwd(old), add = TRUE)
  sys.source(f, envir = env)
  env
}

# Inputs for a small but complete empirical design on pacific_agri, shared by
# the server tests so they differ only in the condition under test.
app_demo_inputs <- function(...) {
  utils::modifyList(list(
    data_source = "builtin",
    builtin_ds  = "Pacific Agriculture (8000 farms)",
    x_vars = c("land_area", "labour_units"),
    y_vars = c("crop_yield", "livestock_value", "expenditure"),
    H_1 = 2, H_2 = 2,
    omega_1 = 1, omega_2 = 1, omega_3 = 1,
    n_total = 120, n_min = 2,
    wmode = "empirical", use_cv = FALSE,
    constraint = "n_total", cost_total = 1500,
    cost_fixed = 0, cost_unit = 1,
    method = "cobyla", n_starts = 1, max_iter = 5,
    run_btn = 1), list(...))
}
