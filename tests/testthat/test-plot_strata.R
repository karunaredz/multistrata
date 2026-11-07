# plot_strata() dispatches on the number of stratification variables:
#   p == 1 -> .plot_strata_1d()   marginal density with boundary cuts
#   p == 2 -> .plot_strata_2d()   boundary map plus a 3-D perspective view
#   p >= 3 -> .plot_strata_nd()   p x p pairs matrix
# Each renderer is exercised on a null graphics device: these are smoke tests
# that the drawing code runs to completion on a real fitted design, not
# comparisons of pixels.

plot_fixture <- function(p, H, N = 900, n_total = 120) {
  set.seed(515)
  pop <- sim_population(N = N, p = p, q = 2, scenario = "independent")
  X <- as.matrix(pop[, paste0("X", seq_len(p)), drop = FALSE])
  Y <- as.matrix(pop[, c("Y1", "Y2")])
  res <- stratify_multi(X, Y, H_vec = H, n_total = n_total,
                        omega = c(0.5, 0.5), n_starts = 2, max_iter = 25,
                        verbose = FALSE)
  list(X = X, Y = Y, res = res)
}

with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
  }, add = TRUE)
  force(code)
}

test_that("plot_strata renders the 1-D case", {
  f <- plot_fixture(p = 1, H = 3)
  with_null_device(
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y)))
})

test_that("plot_strata renders the 2-D case with and without study variables", {
  f <- plot_fixture(p = 2, H = c(2, 3))
  with_null_device({
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y))
    # study_vars is optional for p == 2 (the density surface does not need it)
    expect_silent(plot_strata(f$res, data = f$X))
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y,
                              show_allocation = FALSE))
  })
})

test_that("plot_strata renders the n-D pairs matrix", {
  f <- plot_fixture(p = 3, H = c(2, 2, 2))
  with_null_device(
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y)))
})

test_that("plot_strata honours jitter_n and a custom palette", {
  f <- plot_fixture(p = 2, H = c(2, 2))
  with_null_device({
    # jitter_n below N forces the subsampling branch
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y,
                              jitter_n = 200L))
    expect_silent(plot_strata(f$res, data = f$X, study_vars = f$Y,
                              palette = c("#1B9E77", "#D95F02",
                                          "#7570B3", "#E7298A")))
  })
})

test_that("the plot and print methods run on a fitted design", {
  f <- plot_fixture(p = 2, H = c(2, 2))
  with_null_device({
    expect_silent(plot(f$res, which = "allocation"))
    expect_silent(plot(f$res, which = "weights"))
    expect_error(plot(f$res, which = "nonsense"),
                 "must be 'allocation' or 'weights'")
  })
  expect_output(print(f$res))
  expect_output(print(summary(f$res)))
})

test_that("a rejected 'which' opens no graphics device", {
  f <- plot_fixture(p = 2, H = c(2, 2))

  # par() opens the default device when none is active, which in a
  # non-interactive session writes Rplots.pdf into the working directory.
  # An invalid argument must therefore be rejected before any graphics call.
  while (grDevices::dev.cur() > 1L) grDevices::dev.off()
  before <- grDevices::dev.cur()
  expect_error(plot(f$res, which = "nonsense"),
               "must be 'allocation' or 'weights'")
  expect_identical(grDevices::dev.cur(), before)
})

test_that("plot_strata_rgl runs when rgl is available", {
  skip_on_cran()
  # skip_if_not_installed() calls requireNamespace(), which triggers rgl's
  # .onLoad. Where rgl is installed but cannot load -- recent macOS ships no
  # OpenGL, so the dlopen of libGLU fails -- that emits warnings on its way
  # to returning FALSE, and they surface in the test report and in R CMD
  # check. Probe quietly instead, and say which of the two cases it is.
  rgl_ok <- suppressWarnings(suppressMessages(
    requireNamespace("rgl", quietly = TRUE)))
  if (!rgl_ok)
    skip(if (nzchar(system.file(package = "rgl")))
           "rgl is installed but cannot load (no OpenGL on this platform)"
         else "rgl is not installed")

  f <- plot_fixture(p = 2, H = c(2, 2))
  # rgl needs an interactive OpenGL device; fall back to the null device so
  # this is safe in a headless check.
  withr_null <- options(rgl.useNULL = TRUE)
  on.exit(options(withr_null), add = TRUE)
  # rgl's NULL device cannot set a window title, so it warns about the
  # windowTitle argument. That is a property of the headless device, not of
  # plot_strata_rgl(), so only genuine errors are of interest here.
  expect_no_error(suppressWarnings(plot_strata_rgl(f$res, data = f$X)))
  try(rgl::close3d(), silent = TRUE)
})

test_that("plot_strata_rgl distinguishes rgl absent from rgl unloadable", {
  f <- plot_fixture(p = 2, H = c(2, 2))

  # Case 1: rgl not installed at all -- tell the user to install it.
  testthat::local_mocked_bindings(
    .rgl_available = function() FALSE,
    .rgl_installed = function() FALSE)
  expect_error(plot_strata_rgl(f$res, data = f$X),
               "Package 'rgl' is required", fixed = TRUE)
  expect_error(plot_strata_rgl(f$res, data = f$X),
               'install.packages("rgl")', fixed = TRUE)
})

test_that("plot_strata_rgl explains an unloadable rgl instead of blaming install", {
  f <- plot_fixture(p = 2, H = c(2, 2))

  # Case 2: rgl is installed but will not load. This is the situation on
  # recent macOS, which ships no OpenGL, so rgl's dlopen of libGLU fails.
  # Telling that user to install rgl is useless -- they have it -- so the
  # message must name the real cause and point at the static alternative.
  testthat::local_mocked_bindings(
    .rgl_available = function() FALSE,
    .rgl_installed = function() TRUE)

  err <- tryCatch(plot_strata_rgl(f$res, data = f$X),
                  error = function(e) conditionMessage(e))
  expect_match(err, "installed but cannot be loaded", fixed = TRUE)
  expect_match(err, "OpenGL", fixed = TRUE)
  expect_match(err, "plot_strata(result, data, study_vars)", fixed = TRUE)
  # It must NOT send them off to reinstall something they already have.
  expect_false(grepl("install.packages", err, fixed = TRUE))
})