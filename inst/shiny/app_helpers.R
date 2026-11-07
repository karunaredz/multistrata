`%||%` <- function(a, b) if (is.null(a)) b else a

.distribution_parameter_names <- c(
  "mean", "sd", "meanlog", "sdlog", "shape", "rate", "scale",
  "shape1", "shape2", "min", "max", "df", "location", "a", "b", "c")

.snapshot_app_inputs <- function(input) {
  value <- function(id) input[[id]]
  source_type <- value("data_source")
  weight_mode <- value("wmode")
  file <- value("file_up")
  file_id <- if (!identical(source_type, "upload") || is.null(file)) NULL else
    unname(as.list(file[1, intersect(c("name", "size", "type"), names(file)), drop = FALSE]))

  # How many per-variable input slots to capture.
  #
  # The UI builds H_j / pf_j / pp_j_* over seq_along(input$x_vars) and
  # omega_k / cv_k over seq_along(input$y_vars), so the snapshot has to span
  # whichever of those is longer (and den_p in density mode). A fixed span of
  # four silently dropped every variable past the fourth: the design still
  # fitted correctly, but the downloadable report and the generated R script
  # described a different design, and edits to those variables did not mark
  # the results stale. The floor of four keeps the snapshot shape unchanged
  # for the common small designs.
  n_counts <- c(4L,
                length(value("x_vars") %||% character()),
                length(value("y_vars") %||% character()),
                suppressWarnings(as.integer(value("den_p") %||% 0L)))
  n_slots <- max(n_counts[is.finite(n_counts)])

  variable_inputs <- lapply(seq_len(n_slots), function(j) list(
    H = value(paste0("H_", j)),
    omega = value(paste0("omega_", j)),
    cv = value(paste0("cv_", j)),
    family = if (identical(weight_mode, "parametric")) value(paste0("pf_", j)) else NULL,
    parameters = if (identical(weight_mode, "parametric"))
      setNames(lapply(.distribution_parameter_names, function(nm)
        value(paste0("pp_", j, "_", nm))), .distribution_parameter_names) else NULL,
    density_family = if (identical(source_type, "density")) value(paste0("den_pf_", j)) else NULL,
    density_range = if (identical(source_type, "density"))
      c(value(paste0("den_lo_", j)), value(paste0("den_hi_", j))) else NULL,
    density_H = if (identical(source_type, "density")) value(paste0("den_H_", j)) else NULL,
    density_omega = if (identical(source_type, "density")) value(paste0("den_om_", j)) else NULL,
    density_parameters = if (identical(source_type, "density"))
      setNames(lapply(.distribution_parameter_names, function(nm)
        value(paste0("den_pp_", j, "_", nm))), .distribution_parameter_names) else NULL
  ))

  list(
    data_source = source_type,
    builtin_ds = if (identical(source_type, "builtin")) value("builtin_ds") else NULL,
    uploaded_file = file_id,
    simulation = if (identical(source_type, "simulate"))
      list(N = value("sim_N"), p = value("sim_p"), q = value("sim_q"),
           scenario = value("sim_scenario"), seed = value("sim_seed")) else NULL,
    density = if (identical(source_type, "density"))
      list(p = value("den_p"), N = value("den_N")) else NULL,
    x_vars = value("x_vars"), y_vars = value("y_vars"),
    n_total = value("n_total"), n_min = value("n_min"),
    weight_mode = if (identical(source_type, "density")) "density" else weight_mode,
    constraint = value("constraint"),
    cost_total = value("cost_total"), cost_fixed = value("cost_fixed"),
    cost_unit = value("cost_unit"), use_cv = value("use_cv"),
    method = value("method"), n_starts = value("n_starts"),
    max_iter = value("max_iter"),
    # The seed determines the design, so changing it must mark the displayed
    # results stale just as any other solver setting does.
    run_seed = value("run_seed"), variables = variable_inputs
  )
}

.validate_distribution_spec <- function(spec, range = NULL, label = "distribution") {
  family <- spec$family %||% ""
  supported <- c("norm", "lnorm", "gamma", "exp", "weibull", "beta", "unif",
                 "chisq", "t", "cauchy", "logis", "pareto", "triang")
  if (!family %in% supported) stop(label, ": unsupported distribution family.")
  supplied <- unlist(spec[setdiff(names(spec), "family")], use.names = FALSE)
  if (length(supplied) && any(!is.finite(supplied)))
    stop(label, ": distribution parameters must be finite.")
  positive <- switch(family,
    norm = "sd", lnorm = "sdlog", gamma = c("shape", "rate"),
    exp = "rate", weibull = c("shape", "scale"),
    beta = c("shape1", "shape2"), chisq = "df", t = "df",
    cauchy = "scale", logis = "scale", pareto = c("shape", "scale"),
    character())
  if (length(positive)) {
    vals <- unlist(spec[positive], use.names = FALSE)
    if (length(vals) != length(positive) || any(!is.finite(vals)) || any(vals <= 0))
      stop(label, ": ", paste(positive, collapse = ", "), " must be positive.")
  }
  if (family == "unif" &&
      (!is.finite(spec$min) || !is.finite(spec$max) || spec$max <= spec$min))
    stop(label, ": uniform max must exceed min.")
  if (family == "triang" &&
      (!all(is.finite(c(spec$a, spec$b, spec$c))) || spec$b <= spec$a ||
       spec$c < spec$a || spec$c > spec$b))
    stop(label, ": triangular parameters require a < b and a <= c <= b.")
  if (!is.null(range) &&
      (length(range) != 2L || any(!is.finite(range)) || range[2] <= range[1]))
    stop(label, ": range upper must exceed range lower.")
  invisible(TRUE)
}

.validate_design_setup <- function(N, H_vec, n_total, n_min,
                                   constraint = "n_total", cost_total = NULL,
                                   cost_fixed = 0, cost_unit = 1) {
  if (!is.finite(N) || N < 2) stop("The population must contain at least two units.")
  if (!length(H_vec) || any(!is.finite(H_vec)) || any(H_vec < 2) ||
      any(H_vec != as.integer(H_vec)))
    stop("Each number of strata must be an integer of at least 2.")
  H_total <- prod(H_vec)
  if (H_total > 512L)
    stop("The app limits a single run to 512 total strata; reduce one or more H values.")
  if (H_total > N)
    stop("The requested ", H_total, " strata exceed the ", N, " population units.")
  if (!is.finite(n_min) || n_min < 1 || n_min != as.integer(n_min))
    stop("Minimum sample per stratum must be a positive integer.")
  if (!is.finite(n_total) || n_total < 1 || n_total != as.integer(n_total))
    stop("Total sample size must be a positive integer.")
  minimum_n <- H_total * n_min
  if (constraint == "n_total" && n_total < minimum_n)
    stop("Total sample n must be at least ", minimum_n, " for ", H_total,
         " strata with n_min = ", n_min, ".")
  if (constraint == "n_total" && n_total > N)
    stop("Total sample n cannot exceed the population size N = ", N, ".")
  if (constraint == "cost") {
    if (!all(is.finite(c(cost_total, cost_fixed, cost_unit))) ||
        cost_total <= cost_fixed || cost_unit <= 0)
      stop("Cost mode requires a positive unit cost and a budget above fixed overhead.")
    minimum_cost <- cost_fixed + minimum_n * cost_unit
    if (cost_total < minimum_cost)
      stop("Budget must be at least ", format(minimum_cost),
           " to fund the minimum allocation.")
  }
  invisible(TRUE)
}

.validate_analysis_frame <- function(data, study_vars) {
  if (!is.matrix(data)) data <- as.matrix(data)
  if (!is.matrix(study_vars)) study_vars <- as.matrix(study_vars)
  if (!nrow(data) || nrow(data) != nrow(study_vars))
    stop("Stratification and study data must have the same non-zero row count.")
  if (!is.numeric(data) || !is.numeric(study_vars) ||
      any(!is.finite(data)) || any(!is.finite(study_vars)))
    stop("Selected variables must contain only finite numeric values.")
  constant <- apply(data, 2L, function(x) diff(range(x)) <= sqrt(.Machine$double.eps))
  if (any(constant))
    stop("Stratification variables must vary: ",
         paste(colnames(data)[constant] %||% which(constant), collapse = ", "), ".")
  invisible(TRUE)
}
