test_that("run_app() is exported and validates its port without launching", {
  skip_if_no_app()
  testthat::skip_if_not_installed("jsonlite")
  expect_true(is.function(run_app))
  expect_error(run_app(port = 99999), "1 to 65535")
})

test_that("input snapshots retain only settings relevant to the active source", {
  skip_if_no_app()
  app_env <- load_app_env()
  base <- list(
    data_source = "builtin", builtin_ds = "Pacific Agriculture (8000 farms)",
    file_up = data.frame(name = "ignored.csv", size = 10, type = "text/csv"),
    wmode = "empirical", constraint = "n_total", n_total = 100, n_min = 2,
    use_cv = FALSE, method = "cobyla", n_starts = 1, max_iter = 5)
  snap <- app_env$.snapshot_app_inputs(base)

  expect_identical(snap$data_source, "builtin")
  expect_null(snap$uploaded_file)
  expect_null(snap$simulation)
  expect_null(snap$density)

  changed_hidden <- base; changed_hidden$sim_N <- 9999
  expect_identical(snap, app_env$.snapshot_app_inputs(changed_hidden))

  changed_design <- base; changed_design$n_total <- 101
  expect_false(identical(snap, app_env$.snapshot_app_inputs(changed_design)))

  # Small designs keep the historical four-slot shape.
  expect_length(snap$variables, 4L)
})

test_that("input snapshots span every variable, not just the first four", {
  skip_if_no_app()
  app_env <- load_app_env()

  # Six stratification variables: the snapshot feeds both the staleness check
  # and the downloadable report / generated script, so dropping variables
  # past the fourth would let those disagree with the fitted design.
  wide <- list(
    data_source = "upload", wmode = "empirical", constraint = "n_total",
    n_total = 600, n_min = 2, use_cv = FALSE, method = "cobyla",
    n_starts = 1, max_iter = 5,
    x_vars = paste0("X", 1:6), y_vars = paste0("Y", 1:5),
    H_1 = 2, H_2 = 2, H_3 = 2, H_4 = 2, H_5 = 3, H_6 = 4,
    omega_1 = 1, omega_2 = 1, omega_3 = 1, omega_4 = 1, omega_5 = 2)
  snap <- app_env$.snapshot_app_inputs(wide)

  expect_length(snap$variables, 6L)
  expect_identical(vapply(snap$variables, function(v) v$H %||% NA_real_,
                          numeric(1)),
                   c(2, 2, 2, 2, 3, 4))

  # Editing the fifth variable's H must change the snapshot, which is what
  # makes the app flag stale results.
  bumped <- wide; bumped$H_5 <- 5
  expect_false(identical(snap, app_env$.snapshot_app_inputs(bumped)))

  # ... and so must the fifth study variable's priority weight.
  reweighted <- wide; reweighted$omega_5 <- 3
  expect_false(identical(snap, app_env$.snapshot_app_inputs(reweighted)))
})

test_that("design validation rejects infeasible allocations and budgets", {
  skip_if_no_app()
  app_env <- load_app_env()
  expect_invisible(app_env$.validate_design_setup(1000, c(3, 3), 100, 2))
  expect_error(app_env$.validate_design_setup(1000, c(4, 4), 20, 2),
               "must be at least")
  expect_error(app_env$.validate_design_setup(
    1000, c(3, 3), 100, 2, "cost", cost_total = 20,
    cost_fixed = 5, cost_unit = 1), "Budget must be at least")
  expect_error(app_env$.validate_design_setup(10000, c(8, 8, 8, 2), 3000, 2),
               "512 total strata")
})

test_that("distribution validation catches invalid parameterisations", {
  skip_if_no_app()
  app_env <- load_app_env()
  expect_invisible(app_env$.validate_distribution_spec(
    list(family = "gamma", shape = 2, rate = 1)))
  expect_error(app_env$.validate_distribution_spec(
    list(family = "gamma", shape = -1, rate = 1)), "must be positive")
  expect_error(app_env$.validate_distribution_spec(
    list(family = "triang", a = 0, b = 1, c = 2)), "a <= c <= b")
  expect_error(app_env$.validate_distribution_spec(
    list(family = "norm", mean = 0, sd = 1), c(2, 1)), "range upper")
})

test_that("analysis frame validation rejects non-finite and constant X data", {
  skip_if_no_app()
  app_env <- load_app_env()
  x <- cbind(X1 = 1:10, X2 = seq(2, 20, 2))
  y <- cbind(Y1 = 11:20)
  expect_invisible(app_env$.validate_analysis_frame(x, y))
  expect_error(app_env$.validate_analysis_frame(cbind(X1 = rep(1, 10)), y),
               "must vary")
  x[1, 1] <- Inf
  expect_error(app_env$.validate_analysis_frame(x, y), "finite numeric")
})

test_that("a built-in design runs through the app server", {
  skip_if_no_app()
  app_env <- load_app_env()
  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs())
    session$flushReact()

    expect_null(run_err())
    expect_s3_class(strata_res()$res, "multistrata")
    expect_false(results_stale())
    expect_match(output$code_txt, "data\\(pacific_agri\\)")

    session$setInputs(n_total = 121)
    session$flushReact()
    expect_true(results_stale())
  })
})

test_that("an invalid sample size is reported without replacing a prior result", {
  skip_if_no_app()
  app_env <- load_app_env()
  shiny::testServer(app_env$server, {
    session$setInputs(
      data_source = "builtin",
      builtin_ds = "Business Survey (5000 firms)",
      x_vars = c("employees", "assets"), y_vars = c("turnover", "wage_bill"),
      H_1 = 4, H_2 = 4, omega_1 = 1, omega_2 = 1,
      n_total = 20, n_min = 2, wmode = "empirical", use_cv = FALSE,
      constraint = "n_total", cost_total = 100, cost_fixed = 0, cost_unit = 1,
      method = "cobyla", n_starts = 1, max_iter = 5, run_btn = 1)
    session$flushReact()

    # 16 strata at n_min = 2 needs at least 32 units, so the run must be
    # refused with an explanation rather than producing a broken design.
    expect_null(strata_res())
    expect_match(run_err(), "must be at least 32")
  })
})

test_that("density mode drives benchmark, sample size, strata trace and omega sweep", {
  skip_if_no_app()
  app_env <- load_app_env()
  shiny::testServer(app_env$server, {
    session$setInputs(
      data_source = "density",
      den_p = 2, den_N = 500,
      den_lo_1 = 0, den_hi_1 = 6, den_H_1 = 2, den_om_1 = 0.5,
      den_lo_2 = 0, den_hi_2 = 2, den_H_2 = 2, den_om_2 = 0.5,
      n_total = 80, n_min = 2,
      constraint = "n_total", cost_total = 1000,
      cost_fixed = 0, cost_unit = 1,
      n_starts = 1, max_iter = 5)
    session$flushReact()
    session$setInputs(den_pf_1 = "exp", den_pf_2 = "triang")
    session$flushReact()
    session$setInputs(den_pp_1_rate = 1,
                      den_pp_2_a = 0, den_pp_2_b = 2, den_pp_2_c = 0,
                      run_btn = 1)
    session$flushReact()

    expect_null(run_err())
    expect_true(isTRUE(strata_res()$density))
    expect_s3_class(strata_res()$res, "multistrata")

    b <- bench_res()
    expect_identical(nrow(b), 4L)
    expect_match(paste(b$method, collapse = " "), "density")

    session$setInputs(fn_cv_1 = 20, fn_cv_2 = 20, run_findn = 1)
    session$flushReact()
    expect_s3_class(findn_res()$fn, "findN")
    expect_true(findn_res()$fn$feasible)

    session$setInputs(ch_range = c(2, 3), ch_starts = 1, run_chooseH = 1)
    session$flushReact()
    expect_s3_class(chooseH_res(), "chooseH")
    expect_identical(nrow(chooseH_res()$table), 2L)

    session$setInputs(ox1 = "X1", ox2 = "X2", ogrid = 3, ostarts = 1,
                      run_omega = 1)
    session$flushReact()
    expect_s3_class(omega_res(), "data.frame")
    expect_identical(nrow(omega_res()), 3L)

    expect_match(output$code_txt, "benchmark_density_designs")
    expect_match(output$code_txt, "find_n_density")
    expect_match(output$code_txt, "choose_H_density")
    expect_match(output$code_txt, "omega_sweep_density")
  })
})

test_that("the generated script reproduces the design the app displayed", {
  skip_if_no_app()
  app_env <- load_app_env()
  shown <- NULL; script <- NULL

  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs(
      y_vars = c("crop_yield", "livestock_value"),
      H_1 = 3, H_2 = 2, n_total = 300, n_starts = 4, max_iter = 30,
      run_seed = 2026))
    session$flushReact()
    shown  <<- strata_res()$res
    script <<- output$code_txt
  })

  # The solver draws random initial boundaries, so the script has to carry the
  # seed the app used. Without it the "reproducible R script" tab emits code
  # that produces a different design from the one on screen.
  expect_match(script, "set.seed(2026)", fixed = TRUE)

  # Run the emitted code up to the fit, in a clean environment, and require
  # the design to come back bit-identical.
  src <- strsplit(script, "\n")[[1]]
  upto <- grep("print(res)", src, fixed = TRUE)[1]
  expect_false(is.na(upto))

  sandbox <- new.env(parent = globalenv())
  sandbox$library <- function(...) invisible(NULL)   # already attached
  # The emitted script ends in print(res); swallow it so the suite stays quiet.
  invisible(utils::capture.output(
    eval(parse(text = paste(src[seq_len(upto)], collapse = "\n")),
         envir = sandbox)))

  expect_identical(sandbox$res$objective, shown$objective)
  expect_identical(sandbox$res$boundaries, shown$boundaries)
  expect_identical(sandbox$res$allocation, shown$allocation)
})

test_that("changing the seed marks the results stale", {
  skip_if_no_app()
  app_env <- load_app_env()
  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs(run_seed = 1))
    session$flushReact()
    expect_false(results_stale())
    session$setInputs(run_seed = 7)
    session$flushReact()
    # The seed determines the design, so it must behave like any other
    # solver setting.
    expect_true(results_stale())
  })
})

test_that("changing an input raises the stale-results banner", {
  skip_if_no_app()
  app_env <- load_app_env()
  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs())
    session$flushReact()
    session$setInputs(n_total = 121)
    session$flushReact()
    expect_true(results_stale())
    expect_match(output$ui_stale$html, "Settings changed")
  })
})

test_that("the Deploy tab scores a frame and attaches valid design weights", {
  skip_if_no_app()
  app_env <- load_app_env()
  out <- NULL

  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs(
      y_vars = c("crop_yield", "livestock_value")))
    session$flushReact()
    session$setInputs(dep_src = "self", dep_seed = 1, dep_planned = FALSE,
                      run_deploy = 1)
    session$flushReact()
    out <<- deploy_res()
    expect_false(is.null(output$ui_deploy))
  })

  expect_null(out$err)
  expect_identical(out$n_frame, 8000L)

  # Every planned sample unit is drawn when no stratum is short of frame units.
  expect_identical(out$n_drawn, 120L)
  expect_true(all(out$cover$shortfall == 0L))
  expect_identical(out$cover$n_h_drawn, out$cover$n_h_planned)

  # as_svydesign() must attach the four design columns ...
  expect_true(all(c(".stratum", ".N_h", ".n_h", ".weight") %in%
                    names(out$wframe)))
  # ... with weights exactly N_h / n_h ...
  expect_equal(out$wframe$.weight, out$wframe$.N_h / out$wframe$.n_h,
               tolerance = 1e-12)
  # ... which must sum to the frame population size: the design-unbiasedness
  # check on the Horvitz-Thompson weights.
  expect_equal(sum(out$wframe$.weight), out$n_frame, tolerance = 1e-8)

  expect_identical(nrow(out$scored), out$n_frame)
  expect_true(".stratum" %in% names(out$scored))
})

test_that("the Deploy tab reports a missing or unusable frame instead of failing", {
  skip_if_no_app()
  app_env <- load_app_env()
  err <- NULL

  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs())
    session$flushReact()
    # "upload" selected but no file chosen
    session$setInputs(dep_src = "up", dep_seed = 1, dep_planned = FALSE,
                      run_deploy = 1)
    session$flushReact()
    err <<- deploy_res()$err
  })

  expect_match(err, "Choose a frame file first")
})

test_that("the app degrades rather than fails when 'plotly' is unavailable", {
  skip_if_no_app()

  # plotly is a Suggests: the app must still load, build its server, and
  # produce a design without it. Regression guard for the interactive outputs
  # being *defined* unconditionally: render functions are resolved at
  # definition time, before any req(..., has_plotly) guard inside them runs,
  # so an unconditional definition breaks the entire server rather than just
  # the interactive plots.
  app_env <- load_app_env(without = "plotly")
  expect_false(app_env$has_plotly)
  expect_true(is.function(app_env$server))

  # The inert stand-ins must not evaluate their argument, and must not be
  # mistaken for working plot functions.
  expect_null(app_env$plotly_output("pl_map"))
  expect_error(app_env$plot_ly(), "plotly")

  shiny::testServer(app_env$server, {
    do.call(session$setInputs, app_demo_inputs())
    session$flushReact()

    expect_null(run_err())
    expect_s3_class(strata_res()$res, "multistrata")

    # The Visualise tab falls back to the static package plots instead of
    # going blank, and the rest of the workflow is unaffected.
    expect_false(is.null(output$ui_visualise))
    expect_match(output$code_txt, "data\\(pacific_agri\\)")
  })
})

test_that("plotly verbs resolve to the plotly namespace when it is installed", {
  skip_if_no_app()
  testthat::skip_if_not_installed("plotly")

  # layout() is the one that matters: plotly::layout and graphics::layout are
  # unrelated functions, so a bare layout() call must not fall through to
  # graphics.
  app_env <- load_app_env()
  expect_true(app_env$has_plotly)
  expect_identical(app_env$layout, plotly::layout)
  expect_identical(app_env$plot_ly, plotly::plot_ly)
  expect_identical(app_env$render_plotly, plotly::renderPlotly)
  expect_identical(app_env$plotly_output, plotly::plotlyOutput)
})
