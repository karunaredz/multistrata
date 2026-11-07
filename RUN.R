# =============================================================================
# RUN.R -- guided tour and smoke test for the multistrata package
#
# Exercises every exported function on the shipped example populations, prints
# what each one answers, and finishes by launching the Shiny app.
#
# HOW TO RUN
#   Whole script:      source("RUN.R")          (or: Rscript RUN.R)
#   Section by section: open in RStudio and run the blocks between the ==== rules
#
# The script is read-only: it writes nothing to disk and installs nothing. It
# runs from the installed package if there is one, otherwise from this source
# directory via pkgload. Runtime is roughly 1-2 minutes before the app opens.
#
# A failure is reported with the section it came from and the script stops, so
# the last thing printed is always the thing that broke.
# =============================================================================

# ---- load the package -------------------------------------------------------
if (requireNamespace("multistrata", quietly = TRUE)) {
  library(multistrata)
  message("Using installed multistrata ",
          utils::packageVersion("multistrata"))
} else if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
  message("Loaded multistrata from source: ", normalizePath("."))
} else {
  stop("multistrata is not installed, and this is not the package source ",
       "directory. Either install the package, or setwd() to the folder ",
       "containing DESCRIPTION and re-run.", call. = FALSE)
}

# The solvers start from randomly drawn boundaries, so a seed is what makes
# this run repeatable. Change it to see how stable the designs are.
SEED <- 1
set.seed(SEED)

# Small settings keep the tour quick. Raise n_starts for production designs
# (the package default is 50); at n_starts = 4 expect a little run-to-run
# variation in the last decimals.
NS <- 4
MI <- 40

rule <- function(title)
  cat("\n", strrep("=", 78), "\n== ", title, "\n", strrep("=", 78), "\n", sep = "")
step <- function(...) cat("\n-- ", ..., "\n", sep = "")

ok <- function(label, expr) {
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "error"))
    stop("RUN.R failed at: ", label, "\n  ", conditionMessage(out), call. = FALSE)
  invisible(out)
}


# =============================================================================
# 1. Example populations and simulation  --  sim_population()
# =============================================================================
rule("1. Example populations")

data(pacific_agri)
data(health_survey)
data(business_survey)

step("Three populations ship with the package:")
for (nm in c("pacific_agri", "health_survey", "business_survey")) {
  d <- get(nm)
  cat(sprintf("   %-16s %5d x %d   %s\n", nm, nrow(d), ncol(d),
              paste(names(d), collapse = ", ")))
}

step("sim_population() builds synthetic populations. All seven scenarios:")
for (sc in eval(formals(sim_population)$scenario)) {
  pop <- ok(paste("sim_population, scenario", sc),
            sim_population(N = 500, p = 2, q = 2, scenario = sc, seed = 1))
  cat(sprintf("   %-16s %d x %d\n", sc, nrow(pop), ncol(pop)))
}

# The working population for sections 2-8.
X <- pacific_agri[, c("land_area", "labour_units")]          # stratification
Y <- pacific_agri[, c("crop_yield", "livestock_value", "expenditure")]  # study
OM <- c(0.4, 0.4, 0.2)   # priority weights: production variables matter most


# =============================================================================
# 2. The core design  --  stratify_multi(), print(), summary()
# =============================================================================
rule("2. Joint optimum stratification and allocation")

step("stratify_multi(): boundaries and allocation solved together, 4 x 3 = 12 strata")
set.seed(SEED)
res <- ok("stratify_multi (empirical)",
  stratify_multi(data = X, study_vars = Y, H_vec = c(4, 3), n_total = 600,
                 omega = OM, n_starts = NS, max_iter = MI))
print(res)

step("summary() adds the per-stratum table")
print(summary(res))

step("What the fitted object carries:")
cat("   names(res):", paste(names(res), collapse = ", "), "\n")
cat(sprintf("   objective (Phi*)   %.8f\n", res$objective))
cat(sprintf("   optimality gap     %.4f%%\n", 100 * res$optimality_gap))
cat(sprintf("   converged          %s\n", res$converged))


# =============================================================================
# 3. Parametric weights  --  fit_distributions(), stratify_multi(distributions=)
# =============================================================================
rule("3. Parametric stratum weights")

step("fit_distributions(): pick a marginal family per variable by AIC")
fd <- ok("fit_distributions", fit_distributions(X))
cat("   chosen families:", paste(sprintf("%s = %s", names(attr(fd, "chosen")),
                                         attr(fd, "chosen")), collapse = ", "), "\n")
cat("   top of the AIC table:\n")
print(utils::head(attr(fd, "fit_table")[order(attr(fd, "fit_table")$AIC), ], 4))

step("distributions = \"auto\" runs that fit and uses analytic weights")
set.seed(SEED)
res_auto <- ok("stratify_multi (auto)",
  stratify_multi(X, Y, H_vec = c(4, 3), n_total = 600, omega = OM,
                 n_starts = NS, max_iter = MI, distributions = "auto"))
cat(sprintf("   empirical Phi* = %.8f\n   parametric Phi* = %.8f\n",
            res$objective, res_auto$objective))

step("Or state the families yourself, and use the gradient-based solver")
set.seed(SEED)
res_par <- ok("stratify_multi (slsqp, explicit distributions)",
  stratify_multi(X, Y, H_vec = c(3, 3), n_total = 600, omega = OM,
                 n_starts = NS, max_iter = MI, method = "slsqp",
                 distributions = list(
                   land_area    = list(family = "lnorm", meanlog = 1.2, sdlog = 0.8),
                   labour_units = list(family = "gamma", shape = 2, rate = 0.5))))
cat(sprintf("   SLSQP Phi* = %.8f   (gap %.4f%%)\n",
            res_par$objective, 100 * res_par$optimality_gap))

step("compute_stratum_stats_parametric(): analytic moments at given boundaries")
m <- res_par$.meta
sp <- ok("compute_stratum_stats_parametric",
  compute_stratum_stats_parametric(
    data = as.matrix(X), study_vars = as.matrix(Y),
    B = res_par$boundaries, H_vec = m$H_vec,
    alpha = m$alpha, beta = m$beta, dist_list = m$distributions))
cat("   analytic weights W_h sum to", sprintf("%.6f", sum(sp$W_h)), "\n")


# =============================================================================
# 4. Constraints  --  cost, and coefficient-of-variation targets
# =============================================================================
rule("4. Constraints")

step("Budget instead of a fixed n: cost_total with a fixed overhead")
set.seed(SEED)
res_cost <- ok("stratify_multi (cost constraint)",
  stratify_multi(X, Y, H_vec = c(3, 3), n_total = 600, omega = OM,
                 n_starts = NS, max_iter = MI,
                 constraint = "cost", cost_total = 8000,
                 cost_fixed = 500, cost_per_unit = 12))
cat(sprintf("   n achieved under budget: %d\n", sum(res_cost$allocation)))

step("CV targets: cap the coefficient of variation per study variable")
set.seed(SEED)
res_cv <- ok("stratify_multi (cv_targets)",
  stratify_multi(X, Y, H_vec = c(3, 3), n_total = 900, omega = OM,
                 n_starts = NS, max_iter = MI,
                 cv_targets = c(0.05, 0.05, 0.10)))
cat(sprintf("   n used: %d\n", sum(res_cv$allocation)))


# =============================================================================
# 5. Precision diagnostics
#    design_cv(), composite_variance(), variance_lower_bound(),
#    relative_efficiency(), neyman_alloc()
# =============================================================================
rule("5. Precision diagnostics")

ss <- res$stratum_stats
q  <- ncol(Y)
s2 <- as.matrix(ss[, paste0("sigma2_k", seq_len(q))])
mu <- as.matrix(ss[, paste0("mu_k",     seq_len(q))])

step("design_cv(): achieved CV per study variable")
cv <- ok("design_cv", design_cv(s2, mu, ss$W_h, res$allocation, ss$N_h))
cat(paste(sprintf("   %-16s %.3f%%", names(Y), 100 * as.numeric(cv)),
          collapse = "\n"), "\n")

step("composite_variance(): the objective, recomputed from its parts")
phi <- ok("composite_variance",
          composite_variance(s2, ss$W_h, res$allocation, ss$N_h, OM))
cat(sprintf("   recomputed %.10f  vs reported %.10f  (match: %s)\n",
            phi, res$objective, isTRUE(all.equal(phi, res$objective))))

step("variance_lower_bound(): Cauchy-Schwarz bound, and the gap to it")
lb <- ok("variance_lower_bound",
         variance_lower_bound(s2, ss$W_h, ss$N_h, OM, sum(res$allocation)))
cat(sprintf("   lower bound %.10f   attained %.10f   gap %.4f%%\n",
            lb, res$objective, 100 * (res$objective - lb) / lb))

step("relative_efficiency(): ratio of two composite variances")
cat(sprintf("   coarse 2x2 vs this design: %.3f\n",
            ok("relative_efficiency",
               relative_efficiency(res$objective, res_par$objective))))

step("neyman_alloc(): the allocation sub-problem on its own")
na <- ok("neyman_alloc",
         neyman_alloc(s2, ss$W_h, OM, 600, ss$N_h, n_min = 2))
cat("   n_h:", paste(na, collapse = ", "), " (sum", sum(na), ")\n")


# =============================================================================
# 6. Design questions
#    benchmark_designs(), find_n(), choose_H(), omega_sweep()
# =============================================================================
rule("6. Design questions")

step("benchmark_designs(): what does solving jointly actually buy?")
bm <- ok("benchmark_designs",
  benchmark_designs(X, Y, H_vec = c(4, 3), n_total = 600, omega = OM,
                    result_joint = res))
print(bm)

step("find_n(): smallest n meeting a CV target on every study variable")
fn <- ok("find_n", find_n(X, Y, cv_targets = c(0.05, 0.05, 0.08), result = res))
print(fn)

step("choose_H(): how many strata are worth having?")
ch <- ok("choose_H",
  choose_H(X, Y, n_total = 600, H_each = 2:4, omega = OM,
           n_starts = 2, max_iter = 25, plot = FALSE))
print(ch)

step("omega_sweep(): what does protecting one variable cost another?")
sw <- ok("omega_sweep",
  omega_sweep(X, Y, n_total = 600, result = res,
              pair = c("crop_yield", "livestock_value"),
              grid = c(0.2, 0.5, 0.8), n_starts = 2, max_iter = 25))
print(sw)
cat("   swept:", paste(attr(sw, "pair_names"), collapse = " vs "), "\n")


# =============================================================================
# 7. Deployment  --  assign_strata(), as_svydesign()
# =============================================================================
rule("7. Deployment")

step("assign_strata(): score a frame into the fitted strata")
lab <- ok("assign_strata", assign_strata(res, X))
print(table(lab))

step("Draw a stratified sample at the optimal allocation")
set.seed(SEED)
sid  <- ok("assign_strata (integer ids)", assign_strata(res, X, labels = FALSE))
take <- unlist(lapply(seq_along(res$allocation), function(h) {
  idx <- which(sid == h)
  if (!length(idx)) return(integer(0))
  sample(idx, min(length(idx), as.integer(round(res$allocation[h]))))
}))
cat(sprintf("   drew %d of %d planned units\n", length(take),
            sum(as.integer(round(res$allocation)))))

step("as_svydesign(frame_only = TRUE): design weights attached")
frm <- ok("as_svydesign (frame_only)",
          as_svydesign(res, pacific_agri[take, ], frame_only = TRUE))
print(utils::head(frm[, c(".stratum", ".N_h", ".n_h", ".weight")], 4))
cat(sprintf("   weights sum to %.1f; frame N is %d  (these should agree)\n",
            sum(frm$.weight), nrow(pacific_agri)))

if (requireNamespace("survey", quietly = TRUE)) {
  step("as_svydesign(): a survey design object, ready for estimation")
  des <- ok("as_svydesign", as_svydesign(res, pacific_agri[take, ]))
  print(des)
  cat("   design-based mean of crop_yield:\n")
  print(survey::svymean(~crop_yield, des))
} else {
  step("install.packages('survey') to also get a svydesign object here")
}


# =============================================================================
# 8. Figures  --  plot_strata(), plot(), plot_strata_rgl()
# =============================================================================
rule("8. Figures")

step("plot_strata(): boundary diagram. p=2 gives a map plus a 3-D view")
ok("plot_strata (p = 2)", plot_strata(res, data = X, study_vars = Y))

step("plot(): allocation and weight bar charts")
ok("plot allocation", plot(res, which = "allocation"))
ok("plot weights",    plot(res, which = "weights"))

step("p = 1 and p >= 3 use different renderers")
set.seed(SEED)
pop1 <- sim_population(N = 1500, p = 1, q = 2, scenario = "independent", seed = 1)
r1 <- ok("stratify_multi (p = 1)",
  stratify_multi(pop1[, "X1", drop = FALSE], pop1[, c("Y1", "Y2")],
                 H_vec = 4, n_total = 200, n_starts = 2, max_iter = 25))
ok("plot_strata (p = 1)",
   plot_strata(r1, data = pop1[, "X1", drop = FALSE], study_vars = pop1[, c("Y1","Y2")]))

set.seed(SEED)
pop3 <- sim_population(N = 1500, p = 3, q = 2, scenario = "independent", seed = 1)
r3 <- ok("stratify_multi (p = 3)",
  stratify_multi(pop3[, c("X1","X2","X3")], pop3[, c("Y1","Y2")],
                 H_vec = c(2, 2, 2), n_total = 200, n_starts = 2, max_iter = 25))
ok("plot_strata (p = 3, pairs matrix)",
   plot_strata(r3, data = pop3[, c("X1","X2","X3")], study_vars = pop3[, c("Y1","Y2")]))

if (requireNamespace("rgl", quietly = TRUE) && interactive()) {
  step("plot_strata_rgl(): interactive 3-D view (drag to rotate)")
  ok("plot_strata_rgl", plot_strata_rgl(res, data = X))
} else {
  step("plot_strata_rgl() needs the 'rgl' package and an interactive session")
}


# =============================================================================
# 9. Data-free mode  --  design from the marginal densities alone
#    stratify_multi_density() and its four companions
# =============================================================================
rule("9. Data-free (density) mode")

dists <- list(X1 = list(family = "exp",    rate = 1),
              X2 = list(family = "lnorm",  meanlog = 0, sdlog = 0.5))
rngs  <- list(c(0, 6), c(0, 8))

step("stratify_multi_density(): no data at all, only the margins")
set.seed(SEED)
resd <- ok("stratify_multi_density",
  stratify_multi_density(distributions = dists, H_vec = c(3, 2),
                         n_total = 400, N = 5000, ranges = rngs,
                         omega = c(0.5, 0.5), n_starts = NS, max_iter = MI))
print(resd)

step("The same four questions, answered analytically")
bd <- ok("benchmark_density_designs",
         benchmark_density_designs(result_joint = resd))
print(bd)

fnd <- ok("find_n_density",
          find_n_density(result = resd, cv_targets = c(0.10, 0.10)))
cat(sprintf("   find_n_density: smallest feasible n = %s\n",
            if (isTRUE(fnd$feasible)) fnd$n else "not attainable"))

chd <- ok("choose_H_density",
  choose_H_density(distributions = dists, n_total = 400, N = 5000,
                   H_each = 2:3, ranges = rngs, omega = c(0.5, 0.5),
                   n_starts = 2, max_iter = 25, plot = FALSE))
cat("   choose_H_density recommends:", paste(chd$recommended, collapse = " x "), "\n")

swd <- ok("omega_sweep_density",
  omega_sweep_density(distributions = dists, H_vec = c(3, 2), n_total = 400,
                      N = 5000, ranges = rngs, omega = c(0.5, 0.5),
                      grid = c(0.3, 0.7), n_starts = 2, max_iter = 25))
print(swd)


# =============================================================================
# 10. Coverage check -- did the tour touch every exported function?
# =============================================================================
rule("10. Coverage")

exported <- sort(getNamespaceExports("multistrata"))
exported <- exported[!grepl("^(print|plot|summary)\\.", exported)]
called   <- c("sim_population","stratify_multi","fit_distributions",
              "compute_stratum_stats_parametric","design_cv",
              "composite_variance","variance_lower_bound","relative_efficiency",
              "neyman_alloc","benchmark_designs","find_n","choose_H",
              "omega_sweep","assign_strata","as_svydesign","plot_strata",
              "plot_strata_rgl","stratify_multi_density",
              "benchmark_density_designs","find_n_density","choose_H_density",
              "omega_sweep_density","run_app")
missing <- setdiff(exported, called)
cat(sprintf("   exported functions : %d\n   exercised here     : %d\n",
            length(exported), length(intersect(exported, called))))
if (length(missing)) {
  cat("   NOT exercised      :", paste(missing, collapse = ", "), "\n")
} else {
  cat("   every exported function was exercised.\n")
}

cat("\n", strrep("=", 78), "\n",
    "All sections completed without error.\n",
    strrep("=", 78), "\n", sep = "")


# =============================================================================
# 11. Launch the app
#
# run_app() needs shiny, bslib, DT and jsonlite. plotly adds the interactive
# figures, readxl Excel upload, survey the design objects on the Deploy tab,
# pagedown direct PDF export, and FitVerse an extra fitting engine -- all
# optional, and the app says so where a feature is unavailable.
#
# The app blocks the console until you close it, so it comes last.
# =============================================================================
rule("11. Launching the app")

.ui_needed  <- c("shiny", "bslib", "DT", "jsonlite")
.ui_missing <- .ui_needed[!vapply(.ui_needed, requireNamespace, logical(1),
                                  quietly = TRUE)]

if (length(.ui_missing)) {
  cat("   The app needs:", paste(.ui_missing, collapse = ", "), "\n")
  cat("   Install them, then call run_app():\n\n")
  cat(sprintf("     install.packages(c(%s))\n",
              paste(sprintf('"%s"', .ui_missing), collapse = ", ")))
  cat("     run_app()\n")
} else {
  .opt <- c(plotly = "interactive figures", readxl = "Excel upload",
            survey = "survey design objects on the Deploy tab",
            pagedown = "direct PDF report export",
            FitVerse = "extra distribution-fitting engine")
  .have <- vapply(names(.opt), requireNamespace, logical(1), quietly = TRUE)
  cat("   Optional features:\n")
  for (i in seq_along(.opt))
    cat(sprintf("     [%s] %-9s %s\n", if (.have[i]) "x" else " ",
                names(.opt)[i], .opt[i]))

  cat("\n   Things worth trying in the app:\n")
  cat("     - Optimal Design : load 'Pacific Agriculture', set 4 x 3, run\n")
  cat("     - Visualise      : Boundary Map and the 3-D Density tabs\n")
  cat("     - Benchmark      : the same comparison as section 6 above\n")
  cat("     - Deploy         : score the frame, then download the weighted sample\n")
  cat("     - R Code         : the script there reproduces what you just built\n")
  cat("\n   Opening the app. Close the window (or press Esc here) to return.\n\n")
  run_app()
}

