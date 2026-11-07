# ==============================================================================
# multistrata: Professional Shiny Application
# Simultaneous Multivariate Optimum Stratification and Sample Allocation
#
# Companion app for the multistrata methodology and package.
#
# Implements ONLY what the multistrata package provides:
#   stratify_multi()   : joint MSAO optimisation (empirical / parametric,
#                        COBYLA / SLSQP, n_total / cost / CV constraints)
#   neyman_alloc()     : generalised Neyman (water-filling) allocation
#   design_cv()        : per-variable CVs of a design
#   variance_lower_bound(), composite_variance(), relative_efficiency()
#   benchmark_designs(): joint vs proportional / equal / sequential designs
#   sim_population()   : synthetic multivariate populations
#   plot_strata(), plot(): publication-quality static views
#   datasets: pacific_agri, business_survey, health_survey
#
# Required packages are checked by run_app(). Optional packages add Excel,
# Plotly, FitVerse, and direct PDF-export capabilities.
# ==============================================================================

source("app_helpers.R", local = TRUE)  # bundled beside app.R in inst/shiny

library(shiny)
library(bslib)
library(DT)
library(multistrata)

has_plotly <- requireNamespace("plotly", quietly = TRUE)
has_readxl <- requireNamespace("readxl", quietly = TRUE)
has_mass   <- requireNamespace("MASS",   quietly = TRUE)
has_fitverse <- requireNamespace("FitVerse", quietly = TRUE)
has_pagedown <- requireNamespace("pagedown", quietly = TRUE)  # server-side HTML->PDF

# ------------------------------------------------------------------------------
# Optional 'plotly' integration.
#
# 'plotly' lives in the package's Suggests, so the app must stay usable without
# it: the UI already falls back to the static package plots (see
# output$ui_visualise). The server side needs the same care, because
# the render function is resolved when an output is *defined* -- before the
# req(..., has_plotly) guard inside the expression can ever run. Attaching the
# package with library() would therefore make the whole server error out at
# construction time on a machine without plotly, rather than degrade.
#
# The six plotly verbs the app uses are bound explicitly instead. When plotly
# is present they resolve through its namespace, which also keeps layout()
# pointing at plotly::layout rather than graphics::layout. When it is absent
# they become inert: render_plotly() never forces its argument, so the output
# bodies are not evaluated, and plotly_output() contributes nothing to the UI.
# ------------------------------------------------------------------------------
if (has_plotly) {
  render_plotly <- plotly::renderPlotly
  plotly_output <- plotly::plotlyOutput
  plot_ly       <- plotly::plot_ly
  layout        <- plotly::layout
  add_trace     <- plotly::add_trace
  add_surface   <- plotly::add_surface
} else {
  render_plotly <- function(expr, ...) shiny::renderUI(NULL)  # never forces expr
  plotly_output <- function(...) NULL
  plot_ly <- layout <- add_trace <- add_surface <- function(...)
    stop("The 'plotly' package is required for interactive plots. ",
         "Install it with install.packages(\"plotly\").", call. = FALSE)
}

# ==============================================================================
# Static data
# ==============================================================================

builtin_specs <- list(
  "Pacific Agriculture (8000 farms)" = list(
    ds = "pacific_agri",
    x  = c("land_area", "labour_units"),
    y  = c("crop_yield", "livestock_value", "expenditure"),
    par = list(land_area    = list(family = "lnorm"),
               labour_units = list(family = "gamma"))),
  "Business Survey (5000 firms)" = list(
    ds = "business_survey",
    x  = c("employees", "assets"),
    y  = c("turnover", "wage_bill"),
    par = list(employees = list(family = "lnorm"),
               assets    = list(family = "pareto"))),
  "Health Survey (6000 people)" = list(
    ds = "health_survey",
    x  = c("age", "bmi"),
    y  = c("systolic_bp", "cholesterol"),
    par = list(age = list(family = "beta"),
               bmi = list(family = "norm")))
)

# The 13 families supported by the package's parametric mode
family_choices <- c(
  "Log-Normal"       = "lnorm",   "Gamma"           = "gamma",
  "Normal"           = "norm",    "Exponential"     = "exp",
  "Weibull"          = "weibull", "Beta (rescaled)" = "beta",
  "Uniform"          = "unif",    "Chi-squared"     = "chisq",
  "Student-t"        = "t",       "Cauchy"          = "cauchy",
  "Logistic"         = "logis",   "Pareto"          = "pareto",
  "Triangular"       = "triang"
)

family_params <- list(
  norm    = c(mean = 0, sd = 1),
  lnorm   = c(meanlog = 0, sdlog = 1),
  gamma   = c(shape = 2, rate = 1),
  exp     = c(rate = 1),
  weibull = c(shape = 1.5, scale = 1),
  beta    = c(shape1 = 2, shape2 = 2),
  unif    = c(min = 0, max = 1),
  chisq   = c(df = 4),
  t       = c(df = 8),
  cauchy  = c(location = 0, scale = 1),
  logis   = c(location = 0, scale = 1),
  pareto  = c(shape = 2, scale = 1),
  triang  = c(a = 0, b = 2, c = 0)   # decreasing right-triangular (mode at a)
)

# fitdistr-supported families for the one-click MLE prefill
mle_map <- c(norm = "normal", lnorm = "lognormal", gamma = "gamma",
             exp = "exponential", weibull = "weibull",
             cauchy = "cauchy", logis = "logistic")

pal_strata <- c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC",
                "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
                "#E6AB02","#A6761D","#666666","#8DD3C7","#BEBADA",
                "#FB8072","#80B1D3","#FDB462","#B3DE69")

# ==============================================================================
# CSS: family styling shared with the stratifyR 2.0 / FitVerse apps
# ==============================================================================

app_css <- "
:root {
  --ink:#0F172A; --ink-2:#334155; --muted:#64748B; --faint:#94A3B8;
  --line:#E9EEF5; --line-2:#EEF2F8;
  --bg:#F4FAFF; --card:#FFFFFF;
  --navy:#1E3A5F; --navy-2:#274B78;
  --teal:#12B3A6; --teal-d:#0E9C90; --indigo:#5B6CF0;
  --shadow-sm:0 1px 2px rgba(16,24,40,.05), 0 1px 3px rgba(16,24,40,.06);
  --shadow-md:0 4px 14px rgba(16,24,40,.07);
  --shadow-lg:0 10px 30px rgba(16,24,40,.10);
  --r-sm:9px; --r:13px; --r-lg:18px;
}
html { min-height: 100%; }
body {
  min-height: 100%;
  color: var(--ink);
  background-color: #F4FAFF;
}

/* ── Header ─────────────────────────────────────────────────────────────── */
.navbar, .app-header {
  background-color: rgba(248,252,255,0.90);
  backdrop-filter: blur(10px);
  border: 2px solid rgba(18,179,166,.22);
  border-left: 5px solid rgba(18,179,166,.72);
  padding: 0.28rem 0.85rem; box-shadow: 0 8px 20px rgba(15,23,42,.07);
  display: flex; align-items: center; justify-content: space-between;
  position: sticky; top: 0; z-index: 30;
}
.brand-wrap { display: flex; align-items: center; gap: 12px; }
.brand-logo { height: 44px; width: auto; display: block;
  filter: drop-shadow(0 3px 10px rgba(30,58,95,.20)); }
.brand-name {
  font-size: 1.12rem; font-weight: 800; letter-spacing: -0.5px;
  color: var(--navy); font-style: normal;
  display: inline-flex; align-items: center;
}
.brand-badge {
  background: linear-gradient(135deg, var(--teal), var(--indigo)); color: #fff;
  font-size: 0.56rem; font-weight: 800; padding: 3px 8px; border-radius: 999px;
  margin-left: 9px; vertical-align: middle; letter-spacing: 1px;
}
.brand-sub { color: var(--muted); font-size: 0.66rem; margin-top: 0; }
.header-right { color: var(--faint); font-size: 0.68rem; text-align: right;
  line-height: 1.35; }

/* ── Layout ─────────────────────────────────────────────────────────────── */
.bslib-sidebar-layout > .sidebar {
  background-color: rgba(248,252,255,0.88);
  border-right: 3px solid rgba(30,58,95,.18);
  box-shadow: inset -1px 0 0 rgba(18,179,166,.20), 8px 0 28px rgba(15,23,42,.045);
  backdrop-filter: blur(10px);
}
.bslib-sidebar-layout > .sidebar .sidebar-content { padding: 1.05rem 1.05rem 2rem; }
.bslib-sidebar-layout > .main { background: transparent; padding: 0.16rem 0.85rem 0.65rem; }

.sidebar-section-label {
  position: relative;
  font-size: 0.76rem; font-weight: 950; letter-spacing: 0.72px;
  text-transform: uppercase; color: #fff;
  display: flex; align-items: center; gap: 9px;
  margin: 1.32rem 0 0.72rem 0;
  padding: 0.52rem 0.72rem 0.52rem 0.5rem;
  border-radius: 999px;
  background: linear-gradient(135deg, #1E3A5F 0%, #2563EB 48%, #0E9C90 100%);
  box-shadow: 0 10px 22px rgba(30,58,95,.20), inset 0 1px 0 rgba(255,255,255,.22);
  border: 1px solid rgba(255,255,255,.44);
}
.sidebar-section-label:first-child { margin-top: 0.1rem; }
.sidebar-section-label::after {
  content:''; flex:1; height:1px; background: rgba(255,255,255,.36);
}
.step-badge {
  width: 26px; height: 26px; border-radius: 50%; flex: none;
  background: #fff;
  color: var(--navy); font-size:0.78rem; font-weight:950;
  display:inline-flex; align-items:center; justify-content:center;
  box-shadow: 0 4px 12px rgba(15,23,42,.20);
  border: 2px solid rgba(18,179,166,.34);
}
.sidebar-section-label .help-icon {
  background: rgba(255,255,255,.20); color: #fff; border: 1px solid rgba(255,255,255,.32);
}
.sidebar-section-label .help-icon:hover { background:#fff; color:var(--navy); }

/* ── Buttons ────────────────────────────────────────────────────────────── */
.btn-run {
  background: linear-gradient(135deg, var(--indigo), var(--teal)) !important;
  color: #fff !important; font-weight: 800 !important; font-size: 0.95rem !important;
  border: none !important; width: 100%;
  padding: 0.72rem !important; border-radius: var(--r-sm) !important;
  letter-spacing: 0.3px; margin-top: 1rem;
  box-shadow: 0 6px 18px rgba(91,108,240,.32);
  transition: transform .12s ease, box-shadow .12s ease, filter .12s ease;
}
.btn-run:hover { filter: brightness(1.05); transform: translateY(-1px);
  box-shadow: 0 10px 24px rgba(91,108,240,.42); }
.btn-run:active { transform: translateY(0); }
.btn-mle {
  font-size: 0.68rem !important; padding: 3px 10px !important;
  background: #ECFDF9 !important; color: var(--teal-d) !important;
  border: 1px solid #A7E8DF !important; border-radius: 999px !important;
  font-weight: 700 !important;
}
.btn-mle:hover { background: #DCFBF4 !important; }
.btn-sm, .btn-preset {
  border-radius: var(--r-sm) !important; font-weight: 700 !important;
  font-size: 0.78rem !important;
}
.btn-preset {
  background: #fff !important; color: var(--navy-2) !important;
  border: 1px solid var(--line) !important; padding: 0.5rem 0.7rem !important;
  box-shadow: var(--shadow-sm); width: 100%; text-align: left;
  transition: border-color .12s ease, box-shadow .12s ease;
}
.btn-preset:hover { border-color: var(--teal) !important; box-shadow: var(--shadow-md); }
.preset-row { display: grid; gap: 0.5rem; margin-bottom: 0.3rem; }

/* Friendly data-source choice tabs */
.source-choice-wrap, .source-choice-wrap .shiny-options-group {
  overflow: visible !important;
}
.source-choice-wrap .shiny-options-group { display: grid; gap: 0.44rem; }
.source-choice-wrap .form-check,
.source-choice-wrap .radio {
  position: relative; margin: 0 !important; padding: 0 !important;
  overflow: visible !important;
}
.source-choice-wrap input[type='radio'],
.source-choice-wrap .form-check-input {
  display: none !important;
}
.source-choice-wrap .form-check-label,
.source-choice-wrap .radio label {
  display: block; cursor: pointer; padding: 0 !important; margin: 0 !important;
  overflow: visible !important;
}
.source-choice {
  position: relative; overflow: hidden;
  background: linear-gradient(180deg, rgba(255,255,255,.98), rgba(244,248,255,.96));
  border: 2px solid rgba(192,211,235,.92);
  border-left: 8px solid rgba(91,108,240,.24);
  border-radius: 999px;
  padding: 0.62rem 0.9rem 0.6rem 1.05rem;
  box-shadow: 0 7px 18px rgba(15,23,42,.06);
  transition: transform .12s ease, border-color .12s ease, box-shadow .12s ease, background .12s ease;
}
.source-choice::before {
  content: '';
  position: absolute; inset: 0 auto 0 0; width: 0;
  background: linear-gradient(180deg, #7C3AED, #A855F7);
  transition: width .12s ease;
  border-radius: 999px 0 0 999px;
}
.source-choice-wrap .form-check:hover .source-choice,
.source-choice-wrap .radio:hover .source-choice {
  transform: translateX(2px);
  border-color: rgba(124,58,237,.38);
  box-shadow: 0 12px 24px rgba(76,29,149,.10);
}
.source-choice-wrap .form-check.is-active .source-choice,
.source-choice-wrap .radio.is-active .source-choice,
.source-choice-wrap .form-check:has(input[type='radio']:checked) .source-choice,
.source-choice-wrap .radio:has(input[type='radio']:checked) .source-choice,
.source-choice-wrap .form-check-input:checked + .form-check-label .source-choice {
  border-color: rgba(124,58,237,.82);
  border-left-color: #7C3AED;
  background: linear-gradient(135deg, #4C1D95 0%, #7C3AED 48%, #A855F7 100%);
  box-shadow:
    0 16px 30px rgba(124,58,237,.30),
    0 7px 16px rgba(76,29,149,.18),
    inset 0 1px 0 rgba(255,255,255,.20);
  transform: translateX(6px) scale(1.025);
}
.source-choice-wrap .form-check.is-active .source-choice::before,
.source-choice-wrap .radio.is-active .source-choice::before,
.source-choice-wrap .form-check:has(input[type='radio']:checked) .source-choice::before,
.source-choice-wrap .radio:has(input[type='radio']:checked) .source-choice::before,
.source-choice-wrap .form-check-input:checked + .form-check-label .source-choice::before {
  width: 8px;
}
.source-choice-top {
  display:flex; align-items:center; justify-content:space-between; gap:0.5rem;
}
.source-choice-title { font-size:0.86rem; font-weight:950; color:var(--navy); }
.source-choice-title { font-size:0.94rem; font-weight:1000; letter-spacing:-0.2px; }
.source-choice-sub { display:block; font-size:0.72rem; color:var(--ink-2); line-height:1.32; margin-top:0.14rem; }
.source-choice-wrap .form-check.is-active .source-choice-title,
.source-choice-wrap .form-check.is-active .source-choice-sub,
.source-choice-wrap .radio.is-active .source-choice-title,
.source-choice-wrap .radio.is-active .source-choice-sub,
.source-choice-wrap .form-check:has(input[type='radio']:checked) .source-choice-title,
.source-choice-wrap .form-check:has(input[type='radio']:checked) .source-choice-sub,
.source-choice-wrap .radio:has(input[type='radio']:checked) .source-choice-title,
.source-choice-wrap .radio:has(input[type='radio']:checked) .source-choice-sub,
.source-choice-wrap .form-check-input:checked + .form-check-label .source-choice-title,
.source-choice-wrap .form-check-input:checked + .form-check-label .source-choice-sub {
  color: #fff;
}

/* Beginner/advanced mode, presets and run checklist */
.mode-switch .shiny-options-group {
  display: grid; grid-template-columns: 1fr 1fr; gap: 0.4rem; margin-bottom: 0.45rem;
}
.mode-switch .form-check { padding: 0 !important; margin: 0 !important; }
.mode-switch input[type='radio'] { display: none !important; }
.mode-switch .form-check-label {
  width: 100%; text-align: center; cursor: pointer;
  border: 1px solid rgba(207,224,243,.95); border-radius: 999px;
  padding: 0.42rem 0.55rem; background: #fff;
  color: var(--navy-2); font-weight: 900; font-size: 0.74rem;
  box-shadow: 0 4px 10px rgba(15,23,42,.05);
}
.mode-switch input[type='radio']:checked + .form-check-label {
  background: linear-gradient(135deg, #1E3A5F, #0E9C90);
  color: #fff; border-color: transparent;
  box-shadow: 0 9px 18px rgba(18,179,166,.18);
}
.preset-grid {
  display:grid; grid-template-columns: 1fr 1fr; gap:0.42rem; margin:0.45rem 0 0.55rem;
}
.preset-grid .btn {
  margin:0 !important; width:100%; border-radius:999px !important;
  font-size:0.7rem !important; font-weight:900 !important;
  padding:0.42rem 0.5rem !important;
  background:linear-gradient(180deg,#FFFFFF,#F0F7FF) !important;
  color:var(--navy-2) !important; border:1px solid rgba(207,224,243,.95) !important;
  box-shadow:0 4px 10px rgba(15,23,42,.05);
}
.preset-grid .btn:hover {
  border-color:rgba(18,179,166,.55) !important;
  transform:translateY(-1px);
}
.checklist-card {
  background: rgba(255,255,255,.9); border: 1px solid rgba(18,179,166,.20);
  border-radius: var(--r); padding: 0.62rem 0.72rem; margin: 0.62rem 0 0.45rem;
  box-shadow: 0 6px 16px rgba(15,23,42,.045);
}
.checklist-title { font-size:0.72rem; font-weight:950; color:var(--navy); margin-bottom:0.34rem; }
.check-item { display:flex; align-items:center; gap:0.42rem; font-size:0.72rem; color:var(--ink-2); margin:0.18rem 0; }
.check-dot {
  width:17px; height:17px; border-radius:50%; flex:none;
  display:inline-flex; align-items:center; justify-content:center;
  font-size:0.64rem; font-weight:950;
}
.check-ok .check-dot { background:#DCFCE7; color:#067647; }
.check-wait .check-dot { background:#FFF7D6; color:#A15C07; }
.result-insight {
  background:rgba(236,253,249,.78); border:1px solid rgba(18,179,166,.20);
  border-left:4px solid var(--teal); border-radius:var(--r-sm);
  padding:0.52rem 0.72rem; margin:0 0 0.62rem;
  color:var(--navy-2); font-size:0.8rem; font-weight:750;
}

/* ── Form controls ──────────────────────────────────────────────────────── */
.form-control, .form-select, .selectize-input {
  border-radius: var(--r-sm) !important; border-color: var(--line) !important;
  font-size: 0.86rem !important; box-shadow: none !important;
}
.form-control:focus, .form-select:focus, .selectize-input.focus {
  border-color: var(--teal) !important;
  box-shadow: 0 0 0 3px rgba(18,179,166,.15) !important;
}
.control-label, .shiny-input-container > label, label.control-label {
  font-size: 0.78rem !important; font-weight: 700 !important; color: var(--ink-2) !important;
}
.form-check-label { font-size: 0.82rem; color: var(--ink-2); }
.irs--shiny .irs-bar { background: var(--teal); border-color: var(--teal); }
.irs--shiny .irs-handle > i:first-child { background: var(--teal-d); }

/* ── Metric cards ───────────────────────────────────────────────────────── */
.metrics-row {
  display: grid; grid-template-columns: repeat(6, minmax(95px, 1fr));
  gap: 0.42rem; margin: -0.1rem 0 -0.9rem;
}
.metric-card {
  background-color: rgba(255,255,255,0.88);
  border: 1px solid rgba(18,179,166,.24); border-radius: var(--r-sm);
  padding: 0.38rem 0.62rem 0.34rem; text-align: left; position: relative; overflow: hidden;
  box-shadow: 0 8px 22px rgba(15,23,42,.075); transition: transform .12s ease, box-shadow .12s ease;
  min-height: 50px;
}
.metric-card:hover { transform: translateY(-1px); box-shadow: var(--shadow-md); }
.metric-card::before {
  content:''; position:absolute; left:0; top:0; bottom:0; width:4px;
  background: linear-gradient(var(--teal), var(--indigo));
}
.metric-value { font-size: 0.98rem; font-weight: 800; color: var(--ink);
  letter-spacing: -0.5px; font-variant-numeric: tabular-nums; }
.metric-label { font-size: 0.58rem; font-weight: 700; letter-spacing: 0.55px;
  text-transform: uppercase; color: var(--faint); margin-top: 0; }

/* ── Tabs ───────────────────────────────────────────────────────────────── */
.card { border: 2px solid rgba(18,179,166,.26) !important;
  border-top: 3px solid rgba(18,179,166,.46) !important;
  border-radius: var(--r-lg) !important;
  background-color: rgba(248,252,255,0.84) !important;
  box-shadow: 0 18px 44px rgba(15,23,42,.095), 0 0 0 1px rgba(91,108,240,.08) !important;
  overflow: hidden; margin-top: 0 !important;
  backdrop-filter: blur(8px); }
.bslib-card, .card.bslib-card,
.bslib-card[data-full-screen],
.bslib-card[data-bslib-card-init] {
  margin-top: -1.55rem !important;
}
.card-body, .bslib-card .card-body { padding: 0.46rem 0.58rem 0.62rem !important; }
.nav-tabs, .nav.nav-underline {
  border-bottom: 1px solid rgba(18,179,166,.22) !important;
  padding: 0.42rem 0.48rem 0.36rem; gap: 0.36rem;
  background:
    linear-gradient(135deg, rgba(236,253,249,.92), rgba(239,246,255,.94) 45%, rgba(245,243,255,.92));
  align-items:center;
}
.nav-tabs .nav-link, .nav-underline .nav-link {
  position: relative;
  border: 2px solid rgba(255,255,255,.78) !important;
  color: #163B63 !important;
  background:
    linear-gradient(180deg, rgba(255,255,255,.96), rgba(229,245,255,.94)) !important;
  font-size: 0.84rem; font-weight: 950; letter-spacing: 0.36px;
  text-transform: uppercase;
  padding: 0.48rem 0.76rem !important;
  border-radius: 999px !important; transition: all .12s ease;
  box-shadow:
    0 7px 16px rgba(30,58,95,.10),
    inset 0 -1px 0 rgba(18,179,166,.13);
  outline: 1px solid rgba(207,224,243,.9);
}
.nav-tabs .nav-link:hover {
  color: #072B4D !important;
  background:
    linear-gradient(135deg, #E4FFFA, #EAF2FF) !important;
  border-color: rgba(18,179,166,.58) !important;
  transform: translateY(-2px);
  box-shadow: 0 11px 24px rgba(18,179,166,.20);
}
.nav-tabs .nav-link.active, .nav-underline .nav-link.active {
  color: #fff !important;
  background:
    linear-gradient(135deg, #172554 0%, #2563EB 45%, #0E9C90 100%) !important;
  border-color: transparent !important;
  outline: 2px solid rgba(18,179,166,.26);
  box-shadow:
    0 13px 28px rgba(37,99,235,.28),
    0 4px 12px rgba(18,179,166,.22);
}
.nav-tabs .nav-link.active::after, .nav-underline .nav-link.active::after {
  content: '';
  position: absolute; left: 20%; right: 20%; bottom: -0.32rem;
  height: 4px; border-radius: 999px;
  background: linear-gradient(90deg, #5B6CF0, #12B3A6);
  box-shadow: 0 3px 8px rgba(18,179,166,.30);
}
.nav-pills .nav-link {
  font-size: 0.8rem; font-weight: 850; color: #355372;
  border-radius: 999px; padding: 0.42rem 0.82rem;
  background: linear-gradient(180deg, #FFFFFF, #F0F7FF);
  border: 1px solid rgba(207,224,243,.95);
  box-shadow: 0 4px 10px rgba(15,23,42,.06);
  transition: all .12s ease;
}
.nav-pills .nav-link:hover {
  color: var(--navy);
  border-color: rgba(18,179,166,.46);
  transform: translateY(-1px);
}
.nav-pills .nav-link.active {
  background: linear-gradient(135deg, #1E3A5F, #0E9C90) !important;
  color:#fff !important;
  border-color: transparent;
  box-shadow: 0 8px 18px rgba(30,58,95,.20);
}
.tab-content { padding-top: 0.16rem; }
.tab-pane .row { --bs-gutter-y: 0.2rem; }
.tab-pane .shiny-input-container { margin-bottom: 0.25rem; }
.tab-pane .checkbox { margin-top: 0; margin-bottom: 0.15rem; }

/* Results-first visualisation layout */
.visualise-panel .nav-pills { gap: 0.22rem; margin-bottom: 0.22rem !important; }
.visualise-panel .mini-note { margin-bottom: 0.35rem; }
.vis-controls {
  display:flex; flex-wrap:wrap; align-items:end; gap:0.55rem 0.8rem;
  margin: 0 0 0.25rem;
}
.vis-controls .shiny-input-container { margin-bottom: 0 !important; }
.vis-controls .form-check { margin-bottom: 0.25rem; }
.vis-plot { margin-top: 0.1rem; }
@media (max-width: 1250px) {
  .metrics-row { grid-template-columns: repeat(3, minmax(110px, 1fr)); }
}

/* ── Tables ─────────────────────────────────────────────────────────────── */
table.dataTable thead th {
  background: #F4F7FC; color: var(--navy-2); font-size: 0.78rem; font-weight: 700;
  text-transform: none; letter-spacing: 0.2px; border-bottom: 2px solid var(--line) !important;
}
table.dataTable tbody td { font-size: 0.82rem; color: var(--ink-2);
  font-variant-numeric: tabular-nums; }
table.dataTable tbody tr:hover { background: #F5FAFF !important; }
.dataTables_wrapper .dataTables_paginate .paginate_button.current {
  background: var(--teal) !important; border-color: var(--teal) !important;
  color:#fff !important; border-radius:6px; }

/* ── Alert / state boxes ────────────────────────────────────────────────── */
.err-box {
  background: #FEF3F2; border: 1px solid #FDD9D5; border-left: 4px solid #F04438;
  color: #912018; border-radius: var(--r-sm); padding: 0.85rem 1.05rem;
  font-size: 0.85rem; margin-bottom: 1rem; white-space: pre-wrap; box-shadow: var(--shadow-sm);
}
.ok-box {
  background: #ECFDF3; border: 1px solid #C4EED4; border-left: 4px solid #12B76A;
  color: #05603A; border-radius: var(--r-sm); padding: 0.6rem 0.95rem;
  font-size: 0.82rem; margin-bottom: 0.85rem;
}
.warn-box {
  background: #FFFAEB; border: 1px solid #FEEBC0; border-left: 4px solid #F79009;
  color: #93370D; border-radius: var(--r-sm); padding: 0.6rem 0.95rem;
  font-size: 0.82rem; margin-bottom: 0.85rem;
}

.explain-card {
  background: rgba(248,252,255,.92); border: 1px solid rgba(207,224,243,.95);
  border-radius: var(--r); padding: 0.75rem 0.9rem; margin: 0.45rem 0 0.8rem;
}
.explain-title { font-weight:900; color:var(--navy); font-size:0.82rem; margin-bottom:0.18rem; }
.explain-card p, .explain-card li { font-size:0.78rem; color:var(--ink-2); line-height:1.5; }
.explain-card ul { margin-bottom:0; padding-left:1.05rem; }

.bounds-card {
  background: linear-gradient(180deg,#FFFFFF, #FBFDFF);
  border: 1px solid var(--line); border-radius: var(--r);
  padding: 0.85rem 1.05rem; margin-bottom: 0.7rem;
  font-family: 'SF Mono', Menlo, Consolas, monospace; font-size: 0.85rem;
  color: var(--ink-2); box-shadow: var(--shadow-sm);
}
.bounds-var { font-weight: 800; color: var(--navy-2); font-family: inherit; }

/* ── Intro / cards ──────────────────────────────────────────────────────── */
.intro-section {
  background-color: rgba(255,255,255,0.88);
  border: 1px solid rgba(18,179,166,.22); border-radius: var(--r-lg);
  padding: 1.25rem 1.45rem; margin-bottom: 1.1rem;
  box-shadow: 0 10px 26px rgba(15,23,42,.065);
}
.intro-section-title {
  display: flex; align-items: center; gap: 11px;
  font-size: 1.05rem; font-weight: 800; color: var(--ink); margin-bottom: 0.7rem;
}
.num-badge {
  background: linear-gradient(135deg, var(--indigo), var(--teal)); color: #fff;
  width: 28px; height: 28px; border-radius: 9px;
  display: inline-flex; align-items: center; justify-content: center;
  font-size: 0.82rem; font-weight: 900; flex: none; box-shadow: var(--shadow-sm);
}
.intro-body { font-size: 0.89rem; color: var(--ink-2); line-height: 1.6; }
.intro-formula {
  background: rgba(244,247,252,.86); border: 1px solid rgba(221,232,246,.98); border-radius: var(--r-sm);
  padding: 0.7rem 1rem; font-family: 'SF Mono', Menlo, Consolas, monospace;
  font-size: 0.85rem; margin: 0.65rem 0; color: var(--navy); overflow-x: auto;
}
.intro-ref { font-size: 0.76rem; color: var(--faint); line-height: 1.6;
  padding: 0.8rem 0.4rem; }
.intro-ref-title { font-weight: 900; color: var(--navy-2); margin-bottom: 0.32rem; }
.intro-ref ul { margin: 0.2rem 0 0; padding-left: 1.15rem; }
.intro-ref li { margin-bottom: 0.18rem; }

/* ── Empty state ────────────────────────────────────────────────────────── */
.empty-state { text-align: center; padding: 4rem 1rem; color: var(--muted); }
.empty-icon  { font-size: 3rem; opacity: 0.85; }
.empty-title { font-size: 1.12rem; font-weight: 800; color: var(--ink-2); margin-top: 0.7rem; }
.empty-sub   { font-size: 0.86rem; margin-top: 0.35rem; max-width: 460px;
  margin-left:auto; margin-right:auto; line-height:1.55; }

/* ── Help tooltip ───────────────────────────────────────────────────────── */
.help-icon {
  position: relative;
  display: inline-flex; align-items: center; justify-content: center;
  width: 15px; height: 15px; border-radius: 50%; margin-left: 5px;
  background: var(--line); color: var(--muted); font-size: 0.62rem; font-weight: 800;
  cursor: help; vertical-align: middle; user-select: none;
}
.help-icon:hover { background: var(--teal); color: #fff; }
.help-icon::after {
  content: attr(data-tip);
  display: none;
  position: absolute;
  z-index: 1000;
  left: 50%; top: calc(100% + 10px);
  transform: translateX(-50%);
  min-width: 230px; max-width: 310px;
  white-space: normal;
  text-transform: none; letter-spacing: normal;
  text-align: left; line-height: 1.35;
  padding: 0.58rem 0.7rem;
  border-radius: 12px;
  background: #0F2F4F;
  color: #fff;
  border: 1px solid rgba(255,255,255,.24);
  box-shadow: 0 12px 28px rgba(15,23,42,.22);
  font-size: 0.72rem; font-weight: 700;
}
.help-icon::before {
  content: '';
  display: none;
  position: absolute;
  z-index: 1001;
  left: 50%; top: calc(100% + 4px);
  transform: translateX(-50%);
  border: 6px solid transparent;
  border-bottom-color: #0F2F4F;
}
.help-icon:hover::after,
.help-icon:hover::before { display: block; }

.mini-note { font-size: 0.75rem; color: var(--muted); margin-top: -2px;
  margin-bottom: 8px; line-height: 1.5; }
.dl-row { margin: 0.7rem 0; display: flex; gap: 0.6rem; flex-wrap: wrap; }
.guide-title { font-weight: 800; color: var(--navy); font-size: 0.95rem; }
.guide-body { color: var(--ink-2); font-size: 0.82rem; line-height: 1.55; margin-top: 0.25rem; }
.guide-steps { display:flex; flex-wrap:wrap; gap:0.4rem; }
.guide-pill {
  border: 1px solid var(--line); border-radius: 999px; padding: 0.25rem 0.55rem;
  background: #fff; color: var(--muted); font-size: 0.72rem; font-weight: 700;
}
.guide-pill.active { background: #ECFDF9; border-color: #A7E8DF; color: var(--teal-d); }
.support-card {
  background:rgba(248,250,252,.88); border:1px solid rgba(18,179,166,.18); border-radius:var(--r-sm);
  padding:0.65rem 0.8rem; color:var(--ink-2); font-size:0.78rem;
  line-height:1.5; margin:0.4rem 0 0.65rem;
}
.action-card {
  background-color: rgba(255,255,255,0.88);
  border: 1px solid rgba(18,179,166,.22); border-radius: var(--r);
  padding: 0.8rem 0.95rem; margin-bottom: 0.75rem; box-shadow: var(--shadow-sm);
}
.action-card-title { font-weight: 850; color: var(--navy); margin-bottom: 0.25rem; }
.code-box pre {
  border-radius: var(--r); border: 1px solid var(--line); background: #0F172A !important;
  color: #E2E8F0 !important; padding: 1rem !important; font-size: 0.78rem !important;
  max-height: 68vh; overflow: auto;
}

/* subtle scrollbars */
*::-webkit-scrollbar { width: 10px; height: 10px; }
*::-webkit-scrollbar-thumb { background: #D6DEE9; border-radius: 999px;
  border: 2px solid var(--bg); }
*::-webkit-scrollbar-thumb:hover { background: #B9C4D3; }
"

app_js <- "
function updateSourceChoiceActive() {
  var $wrap = $('.source-choice-wrap');
  if (!$wrap.length) return;
  $wrap.find('.form-check, .radio').removeClass('is-active');
  var $checked = $wrap.find('input[type=\"radio\"]:checked');
  if ($checked.length) $checked.closest('.form-check, .radio').addClass('is-active');
}

$(document).on('click', '.source-choice-wrap .source-choice', function() {
  var $radio = $(this).closest('.form-check, .radio').find('input[type=\"radio\"]');
  if ($radio.length) {
    $radio.prop('checked', true).trigger('change');
    if (window.Shiny) Shiny.setInputValue('data_source', $radio.val(), {priority: 'event'});
    updateSourceChoiceActive();
  }
});

$(document).on('change', '.source-choice-wrap input[type=\"radio\"]', updateSourceChoiceActive);
$(document).on('shiny:connected shiny:value shiny:inputchanged', function() {
  window.setTimeout(updateSourceChoiceActive, 0);
});
$(document).ready(updateSourceChoiceActive);
"

# ==============================================================================
# Helpers
# ==============================================================================

num_cols_of <- function(df)
  names(df)[vapply(df, function(c) is.numeric(c) && sum(is.finite(c)) > 10, TRUE)]

fmt_num <- function(x, d = 4) formatC(x, format = "f", digits = d)

load_builtin <- function(ds_name) {
  e <- new.env(parent = emptyenv())
  utils::data(list = ds_name, package = "multistrata", envir = e)
  get(ds_name, envir = e)
}

empty_state <- function(msg = paste("Set up your design in the sidebar on the",
                                    "left, then press Run, results and",
                                    "interactive charts will appear here."),
                        title = "No results yet",
                        icon = "\U0001F4CA")
  div(class = "empty-state",
      div(class = "empty-icon", icon),
      div(class = "empty-title", title),
      div(class = "empty-sub", msg))

# Small circular "?" that shows a plain-language explanation on hover.
help_icon <- function(text)
  tags$span(class = "help-icon", `aria-label` = text, `data-tip` = text, "?")

# Numbered, step-badged sidebar section heading with an optional hover hint.
section_label <- function(n, title, hint = NULL)
  div(class = "sidebar-section-label",
      tags$span(class = "step-badge", n), title,
      if (!is.null(hint)) help_icon(hint))

source_choice <- function(title, subtitle, help = NULL) {
  div(class = "source-choice",
      div(class = "source-choice-top",
          span(class = "source-choice-title", title)),
      div(class = "source-choice-sub", subtitle))
}

.boundaries_table <- function(s) {
  b <- s$res$boundaries
  do.call(rbind, lapply(seq_along(b), function(j)
    data.frame(variable = s$a$Xv[j], cut = seq_along(b[[j]]),
               boundary = b[[j]], stringsAsFactors = FALSE)))
}

.allocation_table <- function(s) {
  ss <- s$res$stratum_stats
  data.frame(stratum = ss$stratum, N_h = ss$N_h, W_h = ss$W_h,
             n_h = ss$n_h, A_h = ss$A_h, stringsAsFactors = FALSE)
}

.html_escape <- function(x)
  gsub("\"", "&quot;",
       gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(x)))))

.df_to_html <- function(df, digits = 4, max_rows = 80) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return("")
  d <- df[seq_len(min(nrow(df), max_rows)), , drop = FALSE]
  for (nm in names(d))
    if (is.numeric(d[[nm]])) d[[nm]] <- fmt_num(d[[nm]], digits)
  head <- paste(sprintf("<th>%s</th>", .html_escape(names(d))), collapse = "")
  rows <- paste(apply(d, 1, function(rw)
    paste0("<tr>", paste(sprintf("<td>%s</td>", .html_escape(rw)), collapse = ""),
           "</tr>")), collapse = "")
  more <- if (nrow(df) > max_rows)
    sprintf("<p class='note'>Showing first %d of %d rows.</p>", max_rows, nrow(df))
    else ""
  paste0("<table><tr>", head, "</tr>", rows, "</table>", more)
}

# Marginal densities with optimal cut-points (used by the app plot and report).
.plot_density_margins <- function(r) {
  m <- r$.meta; p <- m$p
  build_pdf <- get(".build_pdf", asNamespace("multistrata"))
  op <- par(mfrow = c(p, 1), mar = c(4, 4.2, 2.4, 1)); on.exit(par(op))
  for (j in seq_len(p)) {
    lo <- m$alpha[j]; hi <- m$beta[j]
    spec <- m$distributions[[j]]
    pdf  <- tryCatch(build_pdf(spec, lo, hi), error = function(e) NULL)
    xs <- seq(lo, hi, length.out = 512)
    ys <- if (is.null(pdf)) rep(0, length(xs)) else
            vapply(xs, function(x) tryCatch(pdf(x), error = function(e) 0), 1)
    ys[!is.finite(ys)] <- 0
    cuts <- r$boundaries[[j]]; brks <- c(lo, cuts, hi)
    plot(xs, ys, type = "n", xlab = m$x_names[j], ylab = "density",
         main = sprintf("%s  ~  %s", m$x_names[j], spec$family),
         bty = "l", ylim = c(0, max(ys, 1e-9) * 1.08))
    for (i in seq_len(length(brks) - 1L)) {
      sel <- xs >= brks[i] & xs <= brks[i + 1L]
      if (any(sel))
        polygon(c(brks[i], xs[sel], brks[i + 1L]), c(0, ys[sel], 0),
                col = adjustcolor(pal_strata[(i - 1L) %% length(pal_strata) + 1L],
                                  0.35), border = NA)
    }
    lines(xs, ys, lwd = 2, col = "grey20")
    abline(v = cuts, col = "grey10", lwd = 2)
    for (b in cuts)
      mtext(sprintf("%.3f", b), side = 3, at = b, cex = 0.72, col = "grey10")
  }
}

# Render a base-graphics plot to a self-contained PNG data URI for the report.
.plot_datauri <- function(fn, width = 900, height = 560, res = 108) {
  tf <- tempfile(fileext = ".png")
  grDevices::png(tf, width = width, height = height, res = res)
  tryCatch(fn(), error = function(e) {
    graphics::plot.new()
    graphics::title(main = paste("plot unavailable:", conditionMessage(e)))
  })
  grDevices::dev.off()
  on.exit(unlink(tf), add = TRUE)
  raw <- readBin(tf, "raw", n = file.info(tf)$size)
  paste0("data:image/png;base64,", jsonlite::base64_enc(raw))
}

# Build a self-contained HTML design report from a stored result list.
.build_report_html <- function(s) {
  if (is.null(s)) return("<!doctype html><html><body>No design yet.</body></html>")
  r <- s$res; a <- s$a; m <- r$.meta; st <- s$settings
  ss <- r$stratum_stats
  p <- length(a$Xv); q <- length(a$Yv)

  s2  <- as.matrix(ss[, paste0("sigma2_k", seq_len(q)), drop = FALSE])
  mu  <- as.matrix(ss[, paste0("mu_k",     seq_len(q)), drop = FALSE])
  cvs <- tryCatch(design_cv(s2, mu, ss$W_h, r$allocation, ss$N_h),
                  error = function(e) rep(NA_real_, q))

  esc <- function(x) gsub("<", "&lt;", gsub("&", "&amp;", as.character(x)))
  mode_lbl <- if (identical(m$mode, "density")) "Density (data-free)"
              else if (identical(m$mode, "parametric")) "Parametric" else "Empirical"

  metric <- function(v, l) sprintf(
    "<div class='mc'><div class='mv'>%s</div><div class='ml'>%s</div></div>", v, l)
  metrics <- paste0(
    metric(fmt_num(r$objective, 6), "Composite &Phi;"),
    metric(sprintf("%.2f%%", 100 * r$optimality_gap), "Optimality gap"),
    metric(nrow(ss), "Strata"),
    metric(sum(r$allocation), "Sample n"),
    metric(mode_lbl, "Weight mode"))

  bnd <- paste(vapply(seq_len(p), function(j) sprintf(
    "<tr><td class='v'>%s</td><td class='mono'>[ %s &nbsp;|&nbsp; %s &nbsp;|&nbsp; %s ]</td></tr>",
    esc(a$Xv[j]), fmt_num(m$alpha[j], 3),
    paste(fmt_num(r$boundaries[[j]], 4), collapse = " | "),
    fmt_num(m$beta[j], 3)), ""), collapse = "")

  yhead <- paste(vapply(seq_len(q), function(k)
    sprintf("<th>mean %s</th><th>var %s</th>", esc(a$Yv[k]), esc(a$Yv[k])), ""),
    collapse = "")
  rows <- paste(vapply(seq_len(nrow(ss)), function(h) {
    ycells <- paste(vapply(seq_len(q), function(k) sprintf(
      "<td>%s</td><td>%s</td>", fmt_num(ss[[paste0("mu_k", k)]][h], 3),
      fmt_num(ss[[paste0("sigma2_k", k)]][h], 4)), ""), collapse = "")
    sprintf("<tr><td class='v'>%s</td><td>%d</td><td>%s</td><td><b>%d</b></td>%s<td>%s</td></tr>",
            esc(ss$stratum[h]), ss$N_h[h], fmt_num(ss$W_h[h], 4),
            ss$n_h[h], ycells, fmt_num(ss$A_h[h], 5))
  }, ""), collapse = "")

  cvrow <- paste(vapply(seq_len(q), function(k) sprintf(
    "<tr><td class='v'>%s</td><td>%s</td></tr>", esc(a$Yv[k]),
    if (is.na(cvs[k])) "&mdash;" else sprintf("%.2f%%", 100 * cvs[k])), ""),
    collapse = "")

  analyses <- s$analysis %||% list()
  benchmark_block <- if (is.data.frame(analyses$benchmark)) paste0(
    "<h2>Benchmark</h2>",
    "<p class='note'>Relative efficiency above 1 means the joint optimum is more efficient.",
    if (isTRUE(s$density)) " Density mode uses analytic density benchmarks." else "",
    "</p>",
    .df_to_html(analyses$benchmark, digits = 4)
  ) else ""
  findn_block <- if (!is.null(analyses$findn) && inherits(analyses$findn$fn, "findN")) {
    fn <- analyses$findn$fn
    if (isTRUE(fn$feasible)) paste0(
      "<h2>Sample size search</h2>",
      sprintf("<p class='note'>Smallest feasible sample size: <b>n = %d</b>. CV targets: %s.</p>",
              fn$n, paste(sprintf("%.2f%%", 100 * fn$cv_targets), collapse = ", ")),
      .df_to_html(fn$trace, digits = 4)
    ) else paste0(
      "<h2>Sample size search</h2>",
      "<p class='note'>Targets were not attainable within the search range.</p>",
      .df_to_html(fn$trace, digits = 4)
    )
  } else ""
  chooseH_block <- if (!is.null(analyses$chooseH) && inherits(analyses$chooseH, "chooseH")) {
    ch <- analyses$chooseH
    paste0(
      "<h2>Number of strata trace</h2>",
      sprintf("<p class='note'>Recommended grid: <b>%s</b> (%d total strata).</p>",
              paste(ch$recommended, collapse = " × "), prod(ch$recommended)),
      .df_to_html(ch$table, digits = 4)
    )
  } else ""
  omega_block <- if (is.data.frame(analyses$omega)) paste0(
    "<h2>&omega; Explorer</h2>",
    "<p class='note'>Sensitivity of the design to priority weights for the selected variable pair.</p>",
    .df_to_html(analyses$omega, digits = 4)
  ) else ""

  bench_note <- ""
  if (is.data.frame(analyses$benchmark) &&
      all(c("method", "relative_efficiency") %in% names(analyses$benchmark))) {
    alt <- analyses$benchmark[analyses$benchmark$relative_efficiency > 1, , drop = FALSE]
    if (nrow(alt)) {
      best_i <- which.max(alt$relative_efficiency)
      bench_note <- sprintf(
        "Against the strongest listed alternative, the joint design is %.2fx as efficient (%s).",
        alt$relative_efficiency[best_i], esc(alt$method[best_i]))
    } else {
      bench_note <- "Benchmark alternatives were computed; no simpler alternative improved on the joint design."
    }
  }
  fn_note <- ""
  if (!is.null(analyses$findn) && inherits(analyses$findn$fn, "findN")) {
    fn <- analyses$findn$fn
    fn_note <- if (isTRUE(fn$feasible))
      sprintf("The smallest sample satisfying the selected CV targets is n = %d.", fn$n)
    else
      "The selected CV targets were not attainable within the tested sample-size range."
  }
  ch_note <- ""
  if (!is.null(analyses$chooseH) && inherits(analyses$chooseH, "chooseH")) {
    ch <- analyses$chooseH
    ch_note <- sprintf("The strata trace recommends %s strata (%d total).",
                       paste(ch$recommended, collapse = " \u00d7 "), prod(ch$recommended))
  }
  omega_note <- if (is.data.frame(analyses$omega))
    "\u03c9 Explorer results are included to show sensitivity to priority-weight choices."
  else ""
  summary_items <- c(
    sprintf("This %s design uses %d strata and allocates n = %d observations from N = %s.",
            tolower(mode_lbl), nrow(ss), sum(r$allocation),
            format(m$N, big.mark = ",")),
    sprintf("The achieved composite objective is %s with an optimality gap of %.2f%%.",
            fmt_num(r$objective, 6), 100 * r$optimality_gap),
    bench_note, fn_note, ch_note, omega_note)
  summary_items <- summary_items[nzchar(summary_items)]
  executive_summary <- paste0(
    "<div class='summary-card'><h2>Executive summary</h2><ul>",
    paste(sprintf("<li>%s</li>", summary_items), collapse = ""),
    "</ul></div>")

  # Inline the hex logo (served from www/) so the report is self-contained.
  logo_path <- "www/logo.svg"
  logo <- if (file.exists(logo_path)) {
    tryCatch(paste(readLines(logo_path, warn = FALSE), collapse = "\n"),
             error = function(e) "")
  } else ""
  logo_badge <- if (nzchar(logo)) paste0("<div class='logo'>", logo, "</div>") else ""

  # Clean subtitle: escape only the variable names, keep separators as HTML.
  sub <- paste0(
    p, " stratification variable", if (p > 1) "s" else "", ": ",
    esc(paste(a$Xv, collapse = ", ")),
    "  &middot;  ", q, " study variable", if (q > 1) "s" else "", ": ",
    esc(paste(a$Yv, collapse = ", ")),
    "  &middot;  ", format(Sys.time(), "%Y-%m-%d %H:%M"))

  # Design figures, rendered server-side and embedded as PNG data URIs.
  img1 <- tryCatch(
    if (isTRUE(s$density))
      .plot_datauri(function() .plot_density_margins(r),
                    width = 900, height = 200 * p + 120)
    else
      .plot_datauri(function() plot_strata(r, data = a$dat, study_vars = a$sv),
                    width = 1100, height = 560),
    error = function(e) NA_character_)
  img2 <- tryCatch(
    .plot_datauri(function() plot(r, which = "allocation"), width = 900, height = 440),
    error = function(e) NA_character_)
  figs <- paste0(
    if (!is.na(img1)) paste0("<img class='fig' src='", img1, "' alt='stratification figure'>") else "",
    if (!is.na(img2)) paste0("<img class='fig' src='", img2, "' alt='allocation figure'>") else "")
  figs_block <- if (nzchar(figs)) paste0("<h2>Design figures</h2>", figs) else ""

  paste0(
"<!doctype html><html lang='en'><head><meta charset='utf-8'>",
"<meta name='viewport' content='width=device-width, initial-scale=1'>",
"<title>multistrata: Design report</title><style>",
"body{font-family:Inter,Segoe UI,Helvetica,Arial,sans-serif;color:#0F172A;",
"background:#F5F8FC;margin:0;padding:0 0 3rem}",
".hd{background:linear-gradient(135deg,#1E3A5F,#12B3A6);color:#fff;padding:1.5rem 2rem;",
"display:flex;align-items:center;gap:18px}",
".hd h1{margin:0;font-size:1.5rem;font-weight:800;letter-spacing:-.5px}",
".hd .sub{opacity:.9;font-size:.85rem;margin-top:.3rem}",
".logo{background:#fff;border-radius:14px;padding:7px 9px;flex:none;display:inline-flex;",
"box-shadow:0 4px 14px rgba(0,0,0,.18)}",
".logo svg{height:66px;width:auto;display:block}",
".wrap{max-width:960px;margin:1.4rem auto;padding:0 1.2rem}",
".mrow{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:.8rem;margin-bottom:1.4rem}",
".mc{background:#fff;border:1px solid #E9EEF5;border-left:4px solid #12B3A6;border-radius:12px;padding:.8rem 1rem}",
".mv{font-size:1.3rem;font-weight:800}.ml{font-size:.68rem;text-transform:uppercase;letter-spacing:.6px;color:#64748B;margin-top:2px}",
".summary-card{background:#fff;border:1px solid #DCE8F6;border-radius:16px;padding:1rem 1.15rem;margin:0 0 1.25rem;box-shadow:0 10px 24px rgba(15,23,42,.06)}",
".summary-card h2{margin-top:0}.summary-card ul{margin:.4rem 0 0 1.1rem;padding:0}.summary-card li{margin:.35rem 0;line-height:1.48;font-size:.9rem}",
"h2{font-size:1.02rem;color:#1E3A5F;margin:1.5rem 0 .6rem}",
"table{width:100%;border-collapse:collapse;background:#fff;border:1px solid #E9EEF5;border-radius:12px;overflow:hidden}",
"th,td{padding:.5rem .7rem;font-size:.83rem;text-align:left;border-bottom:1px solid #EEF2F8}",
"th{background:#F4F7FC;color:#274B78;font-weight:700}",
".v{font-weight:700;color:#274B78}.mono{font-family:SFMono-Regular,Consolas,monospace}",
".fig{width:100%;border:1px solid #E9EEF5;border-radius:12px;margin:.3rem 0 1.1rem;background:#fff}",
".note{font-size:.8rem;color:#64748B;line-height:1.55;margin:.2rem 0 .55rem}",
".pdfbtn{position:fixed;top:14px;right:16px;background:#5B6CF0;color:#fff;border:none;",
"padding:.55rem .95rem;border-radius:9px;font-weight:700;font-size:.82rem;cursor:pointer;",
"box-shadow:0 5px 16px rgba(91,108,240,.4);z-index:20}",
".foot{color:#94A3B8;font-size:.75rem;margin-top:1.6rem;text-align:center}",
"@media print{.pdfbtn{display:none}body{background:#fff}",
".hd{-webkit-print-color-adjust:exact;print-color-adjust:exact}",
"table,.fig{page-break-inside:avoid}}",
"</style></head><body>",
"<button class='pdfbtn' onclick='window.print()'>&#8681; Save as PDF</button>",
"<div class='hd'>", logo_badge,
"<div><h1>multistrata &mdash; Optimal design report</h1>",
"<div class='sub'>", sub, "</div></div></div>",
"<div class='wrap'>",
"<div class='mrow'>", metrics, "</div>",
executive_summary,
"<h2>Strata per variable</h2><table><tr><th>Variable</th><th>H</th></tr>",
paste(sprintf("<tr><td class='v'>%s</td><td>%d</td></tr>", esc(a$Xv),
              m$H_vec), collapse = ""), "</table>",
"<h2>Optimal boundaries</h2><table><tr><th>Variable</th><th>[ lower | cut-points | upper ]</th></tr>",
bnd, "</table>",
figs_block,
"<h2>Stratum-level statistics</h2><table><tr><th>Stratum</th>",
"<th>N<sub>h</sub> (population)</th><th>W<sub>h</sub></th>",
"<th>n<sub>h</sub> (sample)</th>", yhead, "<th>A<sub>h</sub></th></tr>",
rows, "</table>",
"<h2>Precision: coefficient of variation</h2><table><tr><th>Study variable</th><th>CV</th></tr>",
cvrow, "</table>",
benchmark_block, findn_block, chooseH_block, omega_block,
"<div class='foot'>Generated by the multistrata companion app &nbsp;&middot;&nbsp; ",
"Total sample n = ", sum(r$allocation),
", population N = ", format(m$N, big.mark = ","), ".</div>",
"</div></body></html>")
}

# Full stratum table: boundary ranges per variable + weights, allocation,
# study-variable moments and A_h. Used by both the on-screen table and the CSV.
.stratum_table <- function(s) {
  r <- s$res; a <- s$a; m <- r$.meta
  ss <- r$stratum_stats
  p <- length(a$Xv); q <- length(a$Yv)
  full_B <- lapply(seq_len(p), function(j) c(m$alpha[j], r$boundaries[[j]], m$beta[j]))
  combos <- expand.grid(lapply(m$H_vec, seq_len), KEEP.OUT.ATTRS = FALSE)

  df <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                   Stratum = ss$stratum)
  for (j in seq_len(p)) {
    idx <- combos[[j]]
    df[[paste0(a$Xv[j], " range")]] <-
      sprintf("[%s, %s]", fmt_num(full_B[[j]][idx], 3),
              fmt_num(full_B[[j]][idx + 1L], 3))
  }
  df[["N_h (population)"]] <- ss$N_h
  df[["W_h (weight)"]]     <- ss$W_h
  df[["n_h (sample)"]]     <- ss$n_h
  for (k in seq_len(q)) {
    df[[paste0("mean ", a$Yv[k])]] <- ss[[paste0("mu_k", k)]]
    df[[paste0("var ",  a$Yv[k])]] <- ss[[paste0("sigma2_k", k)]]
  }
  df[["A_h"]] <- ss$A_h
  df
}

# stratum grid rectangles for a 2-D pair of stratification variables
cell_rects <- function(bx, by, rx, ry) {
  ex <- c(rx[1], bx, rx[2]); ey <- c(ry[1], by, ry[2])
  out <- list(); k <- 0
  for (i in seq_len(length(ex) - 1)) for (j in seq_len(length(ey) - 1)) {
    k <- k + 1
    out[[k]] <- list(x0 = ex[i], x1 = ex[i + 1], y0 = ey[j], y1 = ey[j + 1],
                     cx = (ex[i] + ex[i + 1]) / 2, cy = (ey[j] + ey[j + 1]) / 2,
                     ix = i, iy = j)
  }
  out
}

# ==============================================================================
# UI
# ==============================================================================

ui <- page_sidebar(
  title = NULL,
  theme = bs_theme(version = 5, base_font = font_google("Inter"),
                   heading_font = font_google("Inter")),
  tags$head(
    tags$style(HTML(app_css)),
    tags$script(HTML(app_js))
  ),

  # ── Header ─────────────────────────────────────────────────────────────────
  div(class = "app-header",
    div(class = "brand-wrap",
      tags$img(src = "logo.svg", class = "brand-logo", alt = "multistrata"),
      div(
        div(class = "brand-name", "multistrata",
            tags$span(class = "brand-badge", "MULTI")),
        div(class = "brand-sub",
            "Multivariate Optimum Stratification & Sample Allocation")
      )
    ),
    div(class = "header-right",
        "Research companion app", tags$br(),
        # Read from DESCRIPTION so the header cannot drift out of step with
        # the package version; the literal is only a last resort for the case
        # where the app is opened without the package being loadable.
        paste("multistrata",
              tryCatch(as.character(utils::packageVersion("multistrata")),
                       error = function(e) "1.0-1")))
  ),

  # ── Sidebar ────────────────────────────────────────────────────────────────
  sidebar = sidebar(width = 340, open = "always",

    div(class = "mode-switch",
      radioButtons("ui_mode", NULL,
        choices = c("Beginner" = "beginner", "Advanced" = "advanced"),
        selected = "beginner", inline = FALSE)),

    section_label(1, "Choose your path",
      "Try the app with example population data included, upload your own file, simulate data, or design from distributions only."),
    div(class = "mini-note",
        "Pick the route that matches what you have. The sidebar will then show",
        " only the inputs needed for that route."),
    div(class = "source-choice-wrap",
      radioButtons("data_source", NULL,
        choiceNames = list(
          source_choice("Built-in Dataset",
                        "Try the app with example population data included."),
          source_choice("Upload Data",
                        "Use your own register or population file."),
          source_choice("Simulate Data",
                        "Generate artificial data for testing or teaching."),
          source_choice("Density-based",
                        "Design from assumed distributions, with no data file.")
        ),
        choiceValues = c("builtin", "upload", "simulate", "density"),
        selected = "builtin")),

    conditionalPanel("input.data_source == 'builtin'",
      selectInput("builtin_ds", NULL, choices = names(builtin_specs))),

    conditionalPanel("input.data_source == 'upload'",
      fileInput("file_up", NULL, accept = c(".csv", ".xlsx", ".xls"),
                buttonLabel = "Browse…", placeholder = "CSV or Excel"),
      div(class = "support-card",
          tags$strong("Upload tip: "),
          "Use one row per population unit. Numeric columns can be used as X or Y variables. ",
          "CSV, XLSX and XLS files are accepted; missing or text-only columns will be ignored for design roles."),
      uiOutput("ui_upload_status")),

    conditionalPanel("input.data_source == 'simulate'",
      fluidRow(
        column(6, numericInput("sim_N", "N", 4000, min = 500, max = 50000, step = 500)),
        column(6, selectInput("sim_scenario", "Scenario",
                 c("independent", "correlated", "skewed", "heavy_tail", "mixed"),
                 selected = "skewed"))),
      fluidRow(
        column(4, numericInput("sim_p", "p (X)", 2, min = 1, max = 4)),
        column(4, numericInput("sim_q", "q (Y)", 2, min = 1, max = 4)),
        column(4, numericInput("sim_seed", "Seed", 42, min = 1)))),

    conditionalPanel("input.data_source == 'density'",
      div(class = "mini-note",
          "Build a design directly from marginal distributions, without a",
          " population data file. Use this when you know the approximate",
          " distribution and practical range of each variable."),
      div(class = "support-card",
          tags$strong("Tip: "),
          "Use practical finite lower and upper bounds for each variable. ",
          "These bounds define the population range used by the optimiser."),
      fluidRow(
        column(6, numericInput("den_p", "Variables (X = Y)", 2, min = 1, max = 4)),
        column(6, numericInput("den_N", "Population N", 2000, min = 1, step = 100))),
      uiOutput("ui_density_vars")),

    uiOutput("ui_roles"),
    uiOutput("ui_data_preview"),

    section_label(2, "Design",
      paste("The shape of the design: how many strata (groups) to cut each X",
            "into, the total sample size to spend, and priority weights ω that",
            "say which study variables matter most.")),
    uiOutput("ui_H"),
    fluidRow(
      column(6, numericInput("n_total", "Total sample n", 500, min = 20, step = 10)),
      column(6, numericInput("n_min", "Min per stratum", 2, min = 2, max = 50))),
    conditionalPanel("input.ui_mode == 'advanced'",
      uiOutput("ui_omega")),

    conditionalPanel("input.ui_mode == 'advanced' && input.data_source != 'density'",
      section_label(3, "Stratum weights",
        paste("How the size of each stratum is worked out: Empirical counts the",
              "sample, Parametric uses a fitted distribution per variable, Auto",
              "picks the best-fitting distribution for you.")),
      radioButtons("wmode", NULL,
        choices = c("Empirical (sample proportions)" = "empirical",
                    "Parametric (analytic CDF weights)" = "parametric",
                    "Auto-parametric (AIC best fit)" = "auto"),
        selected = "empirical"),
      conditionalPanel("input.wmode == 'auto'",
        div(class = "mini-note",
            "Families selected per variable by AIC at run time",
            if (has_fitverse) " (FitVerse engine active)." else
            " (builtin engine; install FitVerse to add its engine).")),
      conditionalPanel("input.wmode == 'parametric'",
        uiOutput("ui_parametric"))),
    conditionalPanel("input.ui_mode == 'advanced' && input.data_source == 'density'",
      section_label(3, "Stratum weights",
        paste("How the size of each stratum is worked out: Empirical counts the",
              "sample, Parametric uses a fitted distribution per variable, Auto",
              "picks the best-fitting distribution for you.")),
      div(class = "mini-note",
          "Density mode: weights are analytic (marginal CDF increments) and",
          " stratum moments come from quadrature, set per variable above.")),

    conditionalPanel("input.ui_mode == 'advanced'",
      section_label(4, "Constraint",
        paste("What you hold fixed: a total sample size, or a money budget with a",
              "fixed overhead and a per-unit cost. You can also cap the",
              "coefficient of variation of each study variable.")),
      radioButtons("constraint", NULL,
        choices = c("Fixed sample size (C1)" = "n_total",
                    "Cost / budget (C2)"     = "cost"),
        selected = "n_total"),
      conditionalPanel("input.constraint == 'cost'",
        fluidRow(
          column(6, numericInput("cost_total", "Budget C*", 1500, min = 1)),
          column(6, numericInput("cost_fixed", "Overhead c₀", 0, min = 0))),
        numericInput("cost_unit", "Cost per unit cₕ (scalar)", 1, min = 0.01, step = 0.1)),
      checkboxInput("use_cv", "Add CV targets per study variable (C3)", FALSE),
      conditionalPanel("input.use_cv", uiOutput("ui_cv_targets")),

      div(class = "support-card",
          tags$strong("Advanced solver options: "),
          "Defaults are fine for most uses. More multi-starts and iterations give a more thorough search but take longer."),
      conditionalPanel("input.data_source != 'density'",
        selectInput("method", "Boundary solver",
          c("COBYLA (derivative-free)" = "cobyla",
            "SLSQP (analytic gradient, parametric only)" = "slsqp"))),
      fluidRow(
        column(6, numericInput("n_starts", "Multi-starts", 15, min = 1, max = 100)),
        column(6, numericInput("max_iter", "Max AO iter", 100, min = 5, max = 500))),
      fluidRow(
        column(6, numericInput("run_seed", "Random seed", 1, min = 1, step = 1)),
        column(6, div(class = "mini-note", style = "margin-top:1.7rem;",
          "The solver starts from random boundaries, so the seed fixes the",
          " design. It is written into the generated R script."))),
      div(class = "mini-note",
          "For quick exploration, use fewer starts. For final runs, increase",
          " starts and iterations for a more thorough search.")),

    section_label(5, "Run",
      "Choose a preset if useful, check the setup, then run the optimisation."),
    div(class = "preset-grid",
      actionButton("preset_quick", "Quick run"),
      actionButton("preset_final", "Final run"),
      actionButton("preset_teach", "Teaching demo"),
      actionButton("preset_cost", "Cost design")),
    uiOutput("ui_run_checklist"),

    actionButton("run_btn", "▶  Run Joint Optimisation", class = "btn-run"),
    uiOutput("ui_run_status")
  ),

  # ── Main ───────────────────────────────────────────────────────────────────
  uiOutput("ui_metrics"),
  uiOutput("ui_error"),
  uiOutput("ui_stale"),

  navset_card_tab(id = "main_tabs",

    nav_panel(tagList(icon("book-open"), " Introduction"), value = "intro",
      uiOutput("ui_intro")),

    nav_panel(tagList(icon("table-list"), " Optimal Design"), value = "design",
      uiOutput("ui_design")),

    nav_panel(tagList(icon("chart-area"), " Visualise"),
      uiOutput("ui_visualise")),

    nav_panel(tagList(icon("scale-balanced"), " Benchmark"),
      uiOutput("ui_benchmark")),

    nav_panel(tagList(icon("bullseye"), " Precision (CV)"),
      uiOutput("ui_precision")),

    nav_panel(tagList(icon("magnifying-glass-chart"), " Sample size"),
      uiOutput("ui_findn")),

    nav_panel(tagList(icon("layer-group"), " No. of Strata"),
      uiOutput("ui_chooseH")),

    nav_panel(tagList(icon("sliders"), " ω Explorer"),
      uiOutput("ui_omega_explorer")),

    nav_panel(tagList(icon("truck-fast"), " Deploy"), value = "deploy",
      uiOutput("ui_deploy")),

    nav_panel(tagList(icon("file-arrow-down"), " Report"),
      uiOutput("ui_report")),

    nav_panel(tagList(icon("code"), " R Code"),
      uiOutput("ui_code"))
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {

  # ── Data assembly ──────────────────────────────────────────────────────────
  uploaded_df <- reactive({
    req(input$file_up)
    ext <- tolower(tools::file_ext(input$file_up$name))
    df <- tryCatch({
      if (ext == "csv") utils::read.csv(input$file_up$datapath)
      else if (ext %in% c("xlsx", "xls")) {
        if (!has_readxl) stop("Install the 'readxl' package to read Excel files.")
        as.data.frame(readxl::read_excel(input$file_up$datapath))
      } else stop("Unsupported file type: .", ext)
    }, error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    df
  })

  dataset <- reactive({
    switch(input$data_source,
      builtin  = load_builtin(builtin_specs[[input$builtin_ds]]$ds),
      upload   = { req(uploaded_df()); uploaded_df() },
      simulate = sim_population(N = input$sim_N, p = input$sim_p, q = input$sim_q,
                                scenario = input$sim_scenario, seed = input$sim_seed),
      density  = NULL)   # density mode is data-free; roles/H/omega come from ui_density_vars
  })

  output$ui_upload_status <- renderUI({
    if (!identical(input$data_source, "upload")) return(NULL)
    if (is.null(input$file_up)) {
      return(div(class = "mini-note",
        "No file selected yet. Upload a CSV, XLSX or XLS file, then choose the X and Y variables."))
    }
    df <- tryCatch(uploaded_df(), error = function(e) NULL)
    if (is.null(df)) {
      return(div(class = "warn-box",
        "I could not read this file. Please check that it is a CSV, XLSX or XLS file."))
    }
    nc <- num_cols_of(df)
    if (!length(nc)) {
      return(div(class = "warn-box",
        "The file loaded, but I could not find usable numeric columns for the design."))
    }
    div(class = "ok-box",
      sprintf("✓ File loaded: %s rows and %d usable numeric columns.",
              format(nrow(df), big.mark = ","), length(nc)))
  })

  observeEvent(input$preset_quick, {
    updateRadioButtons(session, "ui_mode", selected = "beginner")
    updateRadioButtons(session, "constraint", selected = "n_total")
    updateCheckboxInput(session, "use_cv", value = FALSE)
    updateNumericInput(session, "n_starts", value = 8)
    updateNumericInput(session, "max_iter", value = 80)
    showNotification("Quick run preset applied: faster settings for exploration.",
                     type = "message", duration = 3)
  })

  observeEvent(input$preset_final, {
    updateRadioButtons(session, "ui_mode", selected = "advanced")
    updateRadioButtons(session, "constraint", selected = "n_total")
    updateNumericInput(session, "n_starts", value = 25)
    updateNumericInput(session, "max_iter", value = 200)
    showNotification("Final run preset applied: more thorough optimisation settings.",
                     type = "message", duration = 3)
  })

  observeEvent(input$preset_teach, {
    updateRadioButtons(session, "ui_mode", selected = "beginner")
    updateRadioButtons(session, "data_source", selected = "simulate")
    updateNumericInput(session, "sim_N", value = 2500)
    updateNumericInput(session, "sim_p", value = 2)
    updateNumericInput(session, "sim_q", value = 2)
    updateSelectInput(session, "sim_scenario", selected = "skewed")
    updateNumericInput(session, "n_total", value = 300)
    updateRadioButtons(session, "constraint", selected = "n_total")
    updateNumericInput(session, "n_starts", value = 8)
    updateNumericInput(session, "max_iter", value = 80)
    showNotification("Teaching demo preset applied: simulated skewed population.",
                     type = "message", duration = 3)
  })

  observeEvent(input$preset_cost, {
    updateRadioButtons(session, "ui_mode", selected = "advanced")
    updateRadioButtons(session, "constraint", selected = "cost")
    updateNumericInput(session, "cost_total", value = 5000)
    updateNumericInput(session, "cost_fixed", value = 100)
    updateNumericInput(session, "cost_unit", value = 8)
    showNotification("Cost design preset applied: advanced budget controls are visible.",
                     type = "message", duration = 3)
  })

  output$ui_run_checklist <- renderUI({
    source <- input$data_source %||% "builtin"
    is_density <- identical(source, "density")
    df <- tryCatch(dataset(), error = function(e) NULL)
    data_ok <- is_density || !is.null(df)
    vars_ok <- if (is_density) {
      is.finite(as.numeric(input$den_p)) && as.numeric(input$den_p) >= 1
    } else {
      length(input$x_vars %||% character()) >= 1 &&
        length(input$y_vars %||% character()) >= 1 &&
        length(intersect(input$x_vars %||% character(), input$y_vars %||% character())) == 0
    }
    Hs <- if (is_density) {
      p <- as.integer(input$den_p %||% 0)
      if (p < 1) integer(0) else vapply(seq_len(p), function(j) {
        v <- input[[paste0("den_H_", j)]]
        if (is.null(v) || !is.finite(as.numeric(v))) 3L else as.integer(v)
      }, integer(1))
    } else {
      Xv <- input$x_vars %||% character()
      if (!length(Xv)) integer(0) else vapply(seq_along(Xv), function(j) {
        v <- input[[paste0("H_", j)]]
        if (is.null(v) || !is.finite(as.numeric(v))) 3L else as.integer(v)
      }, integer(1))
    }
    n_total <- as.numeric(input$n_total %||% NA_real_)
    n_min <- as.numeric(input$n_min %||% NA_real_)
    design_ok <- length(Hs) >= 1 && is.finite(n_total) && is.finite(n_min) &&
      n_total >= n_min * prod(Hs)
    ready <- data_ok && vars_ok && design_ok
    item <- function(ok, text) {
      div(class = paste("check-item", if (isTRUE(ok)) "check-ok" else "check-wait"),
          span(class = "check-dot", if (isTRUE(ok)) "✓" else "!"),
          span(text))
    }
    div(class = "checklist-card",
      div(class = "checklist-title", "Setup checklist"),
      item(data_ok, if (is_density) "Density setup selected." else "Population data is ready."),
      item(vars_ok, if (is_density) "Density variables are ready." else "X and Y variables are selected."),
      item(design_ok, "Sample size is feasible for the number of strata."),
      item(ready, if (ready) "Ready to run." else "Fix the highlighted item before running."))
  })

  # ── Variable roles ─────────────────────────────────────────────────────────
  output$ui_roles <- renderUI({
    df <- tryCatch(dataset(), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    nc <- num_cols_of(df)
    def_x <- def_y <- NULL
    if (input$data_source == "builtin") {
      sp <- builtin_specs[[input$builtin_ds]]
      def_x <- sp$x; def_y <- sp$y
    } else if (input$data_source == "simulate") {
      def_x <- grep("^X", nc, value = TRUE); def_y <- grep("^Y", nc, value = TRUE)
    } else {
      def_x <- head(nc, 2); def_y <- head(setdiff(nc, def_x), 2)
    }
    tagList(
      selectizeInput("x_vars", "Stratification variables (X)",
        choices = nc, selected = def_x, multiple = TRUE,
        options = list(maxItems = 4)),
      selectizeInput("y_vars", "Study variables (Y)",
        choices = nc, selected = def_y, multiple = TRUE,
        options = list(maxItems = 4))
    )
  })

  output$ui_data_preview <- renderUI({
    df <- tryCatch(dataset(), error = function(e) NULL)
    req(df, input$x_vars, input$y_vars)
    div(class = "mini-note",
        sprintf("N = %s units · p = %d stratification · q = %d study",
                format(nrow(df), big.mark = ","),
                length(input$x_vars), length(input$y_vars)))
  })

  # ── Dynamic H, omega, CV, parametric inputs ────────────────────────────────
  output$ui_H <- renderUI({
    req(input$x_vars)
    tagList(
      tags$label(class = "control-label", "Strata per variable (Hⱼ)"),
      fluidRow(lapply(seq_along(input$x_vars), function(j)
        column(max(3, floor(12 / length(input$x_vars))),
          numericInput(paste0("H_", j), input$x_vars[j], value = 3,
                       min = 2, max = 8)))),
      uiOutput("ui_H_note")
    )
  })

  output$ui_H_note <- renderUI({
    req(input$x_vars)
    Hs <- vapply(seq_along(input$x_vars),
                 function(j) { v <- input[[paste0("H_", j)]]; if (is.null(v)) 3 else v },
                 numeric(1))
    div(class = "mini-note",
        sprintf("Total strata H = %s", paste(Hs, collapse = " × ")),
        sprintf(" = %d", prod(Hs)))
  })

  output$ui_omega <- renderUI({
    req(input$y_vars)
    q <- length(input$y_vars)
    tagList(
      tags$label(class = "control-label", "Priority weights ω (normalised on run)"),
      fluidRow(lapply(seq_len(q), function(k)
        column(max(3, floor(12 / q)),
          numericInput(paste0("omega_", k), input$y_vars[k],
                       value = round(1 / q, 2), min = 0.01, max = 1, step = 0.05))))
    )
  })

  output$ui_cv_targets <- renderUI({
    req(input$y_vars)
    tagList(
      fluidRow(lapply(seq_along(input$y_vars), function(k)
        column(max(3, floor(12 / length(input$y_vars))),
          numericInput(paste0("cv_", k), paste0("γ ", input$y_vars[k]),
                       value = 0.05, min = 0.001, max = 1, step = 0.005)))),
      div(class = "mini-note",
          "The allocation is adjusted so every selected study variable aims",
          " to stay at or below its CV target.")
    )
  })

  output$ui_parametric <- renderUI({
    req(input$x_vars)
    def_fams <- rep("lnorm", length(input$x_vars))
    if (input$data_source == "builtin") {
      sp <- builtin_specs[[input$builtin_ds]]
      for (j in seq_along(input$x_vars)) {
        pf <- sp$par[[input$x_vars[j]]]
        if (!is.null(pf)) def_fams[j] <- pf$family
      }
    }
    accordion(open = FALSE, !!!lapply(seq_along(input$x_vars), function(j)
      accordion_panel(
        title = paste0("Fⱼ for ", input$x_vars[j]),
        selectInput(paste0("pf_", j), NULL, family_choices,
                    selected = def_fams[j]),
        uiOutput(paste0("ui_pp_", j)),
        if (has_mass) actionButton(paste0("mle_", j),
          "⚙ Estimate by MLE (fitdistr)", class = "btn-mle")
      )))
  })

  # per-variable parameter inputs + MLE prefill
  lapply(1:4, function(j) {
    output[[paste0("ui_pp_", j)]] <- renderUI({
      req(input$x_vars, length(input$x_vars) >= j, input[[paste0("pf_", j)]])
      pd <- family_params[[input[[paste0("pf_", j)]]]]
      fluidRow(lapply(seq_along(pd), function(i)
        column(max(4, floor(12 / length(pd))),
          numericInput(paste0("pp_", j, "_", names(pd)[i]),
                       names(pd)[i], value = unname(pd[i]), step = 0.05))))
    })

    observeEvent(input[[paste0("mle_", j)]], {
      req(input$x_vars, length(input$x_vars) >= j)
      fam <- input[[paste0("pf_", j)]]
      df  <- dataset(); x <- df[[input$x_vars[j]]]; x <- x[is.finite(x)]
      if (fam %in% names(mle_map)) {
        ft <- tryCatch(MASS::fitdistr(x, mle_map[[fam]]), error = function(e) NULL)
        if (is.null(ft)) {
          showNotification("MLE failed for this family/data.", type = "warning")
        } else {
          est <- ft$estimate
          # map fitdistr names to package parameter names
          nm <- names(est)
          nm[nm == "meanlog"] <- "meanlog"; nm[nm == "sdlog"] <- "sdlog"
          for (i in seq_along(est))
            updateNumericInput(session, paste0("pp_", j, "_", nm[i]),
                               value = round(unname(est[i]), 4))
          showNotification(paste0("MLE fitted for ", input$x_vars[j]), type = "message")
        }
      } else if (fam == "unif") {
        updateNumericInput(session, paste0("pp_", j, "_min"), value = round(min(x), 4))
        updateNumericInput(session, paste0("pp_", j, "_max"), value = round(max(x), 4))
      } else if (fam == "triang") {
        updateNumericInput(session, paste0("pp_", j, "_a"), value = round(min(x), 4))
        updateNumericInput(session, paste0("pp_", j, "_b"), value = round(max(x), 4))
        updateNumericInput(session, paste0("pp_", j, "_c"),
                           value = round(as.numeric(names(sort(-table(round(x, 1))))[1]), 4))
      } else if (fam == "pareto") {
        sc <- min(x[x > 0]); sh <- length(x) / sum(log(x[x > 0] / sc))
        updateNumericInput(session, paste0("pp_", j, "_scale"), value = round(sc, 4))
        updateNumericInput(session, paste0("pp_", j, "_shape"), value = round(sh, 4))
      } else {
        showNotification("MLE prefill not available for this family; enter parameters manually.",
                         type = "warning")
      }
    })
  })

  # keep the solver menu honest: SLSQP needs parametric mode
  observe({
    if (identical(input$wmode, "empirical") && identical(input$method, "slsqp")) {
      updateSelectInput(session, "method", selected = "cobyla")
      showNotification("SLSQP requires the parametric mode (analytic weight gradient).",
                       type = "warning")
    }
  })

  # ── Density-mode per-variable inputs (family, params, range, Hⱼ, ωⱼ) ────────
  output$ui_density_vars <- renderUI({
    p <- as.integer(input$den_p %||% 2)
    if (is.na(p) || p < 1) return(NULL)
    accordion(open = TRUE, !!!lapply(seq_len(p), function(j)
      accordion_panel(
        title = paste0("Variable X", j),
        selectInput(paste0("den_pf_", j), "Family", family_choices,
                    selected = if (j == 2) "gamma" else "lnorm"),
        div(class = "mini-note",
            "Choose the marginal distribution and parameter values for this variable."),
        uiOutput(paste0("den_ui_pp_", j)),
        div(class = "support-card",
            tags$strong("Support range: "),
            "the optimiser only uses values inside this range. ",
            "Pick realistic bounds for the population you expect to survey."),
        fluidRow(
          column(6, numericInput(paste0("den_lo_", j), "Range lower",
                                 value = 0, step = 0.5)),
          column(6, numericInput(paste0("den_hi_", j), "Range upper",
                                 value = 10, step = 0.5))),
        fluidRow(
          column(6, numericInput(paste0("den_H_", j), "Strata Hⱼ",
                                 value = 3, min = 2, max = 8)),
          column(6, numericInput(paste0("den_om_", j), "Weight ωⱼ",
                                 value = round(1 / p, 2), min = 0.01,
                                 max = 1, step = 0.05))))))
  })

  # dynamic parameter fields per density variable
  lapply(1:4, function(j) {
    output[[paste0("den_ui_pp_", j)]] <- renderUI({
      req(input$den_p, as.integer(input$den_p) >= j, input[[paste0("den_pf_", j)]])
      pd <- family_params[[input[[paste0("den_pf_", j)]]]]
      fluidRow(lapply(seq_along(pd), function(i)
        column(max(4, floor(12 / length(pd))),
          numericInput(paste0("den_pp_", j, "_", names(pd)[i]),
                       names(pd)[i], value = unname(pd[i]), step = 0.05))))
    })
  })

  # ── Assemble arguments & run ───────────────────────────────────────────────
  gather_args <- function() {
    df <- dataset()
    Xv <- input$x_vars; Yv <- input$y_vars
    validate(need(length(Xv) >= 1, "Select at least one stratification variable."),
             need(length(Yv) >= 1, "Select at least one study variable."),
             need(length(intersect(Xv, Yv)) == 0,
                  "X and Y variable sets must not overlap."))
    dat <- df[, Xv, drop = FALSE]; sv <- df[, Yv, drop = FALSE]
    keep <- stats::complete.cases(dat) & stats::complete.cases(sv)
    dat <- dat[keep, , drop = FALSE]; sv <- sv[keep, , drop = FALSE]

    Hs <- vapply(seq_along(Xv), function(j) {
      v <- input[[paste0("H_", j)]]
      if (is.null(v) || !is.finite(as.numeric(v))) 3L else as.integer(v)
    }, 1L)
    om <- vapply(seq_along(Yv), function(k) {
      v <- input[[paste0("omega_", k)]]
      if (is.null(v) || !is.finite(as.numeric(v))) 1 / length(Yv) else as.numeric(v)
    }, 1)
    om <- om / sum(om)

    dists <- NULL; auto_tbl <- NULL; auto_used <- FALSE
    if (identical(input$wmode, "auto")) {
      if (!"fit_distributions" %in% getNamespaceExports("multistrata"))
        stop("The installed multistrata predates the auto-parametric ",
             "mode. Restart R (Session > Restart R), then re-source ",
             "run_app.R so the updated package is loaded.")
      dists <- multistrata::fit_distributions(dat, verbose = FALSE)
      auto_tbl <- attr(dists, "fit_table")
      attr(dists, "fit_table") <- NULL
      attr(dists, "chosen")    <- NULL
      auto_used <- TRUE
    }
    if (identical(input$wmode, "parametric")) {
      dists <- lapply(seq_along(Xv), function(j) {
        fam <- input[[paste0("pf_", j)]]
        pd  <- family_params[[fam]]
        pars <- lapply(names(pd), function(nm) {
          v <- input[[paste0("pp_", j, "_", nm)]]
          if (is.null(v)) unname(pd[nm]) else as.numeric(v)
        })
        names(pars) <- names(pd)
        c(list(family = fam), pars)
      })
      names(dists) <- Xv
      for (j in seq_along(dists))
        .validate_distribution_spec(dists[[j]], label = paste("Variable", Xv[j]))
    }

    cvt <- NULL
    if (isTRUE(input$use_cv))
      cvt <- vapply(seq_along(Yv), function(k) as.numeric(input[[paste0("cv_", k)]]), 1)

    .validate_analysis_frame(dat, sv)
    .validate_design_setup(
      N = nrow(dat), H_vec = Hs, n_total = as.numeric(input$n_total),
      n_min = as.numeric(input$n_min), constraint = input$constraint %||% "n_total",
      cost_total = input$cost_total, cost_fixed = input$cost_fixed,
      cost_unit = input$cost_unit)

    list(dat = dat, sv = sv, Hs = Hs, om = om, dists = dists, cvt = cvt,
         Xv = Xv, Yv = Yv, auto = auto_used, auto_tbl = auto_tbl)
  }

  # Density (data-free) argument assembly. Study vars = stratification vars.
  gather_density_args <- function() {
    p <- as.integer(input$den_p)
    validate(need(!is.na(p) && p >= 1, "Set at least one variable."))
    Xv <- paste0("X", seq_len(p))

    dists <- lapply(seq_len(p), function(j) {
      fam <- input[[paste0("den_pf_", j)]]
      pd  <- family_params[[fam]]
      pars <- lapply(names(pd), function(nm) {
        v <- input[[paste0("den_pp_", j, "_", nm)]]
        if (is.null(v)) unname(pd[nm]) else as.numeric(v)
      })
      names(pars) <- names(pd)
      c(list(family = fam), pars)
    })
    names(dists) <- Xv

    ranges <- lapply(seq_len(p), function(j)
      c(as.numeric(input[[paste0("den_lo_", j)]]),
        as.numeric(input[[paste0("den_hi_", j)]])))
    Hs <- vapply(seq_len(p), function(j) {
      v <- input[[paste0("den_H_", j)]]
      if (is.null(v) || !is.finite(as.numeric(v))) 3L else as.integer(v)
    }, 1L)
    om <- vapply(seq_len(p), function(j) {
      v <- input[[paste0("den_om_", j)]]
      if (is.null(v) || !is.finite(as.numeric(v))) 1 / p else as.numeric(v)
    }, 1)
    om <- om / sum(om)

    for (j in seq_len(p))
      validate(need(is.finite(ranges[[j]][1]) && is.finite(ranges[[j]][2]) &&
                    ranges[[j]][2] > ranges[[j]][1],
                    sprintf("Variable X%d: range upper must exceed lower.", j)))

    for (j in seq_len(p))
      .validate_distribution_spec(dists[[j]], ranges[[j]], paste("Variable", Xv[j]))
    .validate_design_setup(
      N = as.numeric(input$den_N), H_vec = Hs,
      n_total = as.numeric(input$n_total), n_min = as.numeric(input$n_min),
      constraint = input$constraint %||% "n_total", cost_total = input$cost_total,
      cost_fixed = input$cost_fixed, cost_unit = input$cost_unit)

    list(dists = dists, ranges = ranges, Hs = Hs, om = om,
         N = as.integer(input$den_N), Xv = Xv, Yv = Xv,
         dat = NULL, sv = NULL, auto = FALSE, auto_tbl = NULL, cvt = NULL)
  }

  strata_res <- reactiveVal(NULL)
  run_err    <- reactiveVal(NULL)

  design_snapshot <- reactive(.snapshot_app_inputs(input))
  results_stale <- reactive({
    s <- strata_res()
    !is.null(s) && !identical(s$snapshot, design_snapshot())
  })
  require_current_result <- function() {
    if (!isTRUE(results_stale())) return(TRUE)
    showNotification("Settings changed. Re-run the design before this analysis.",
                     type = "warning", duration = 6)
    FALSE
  }

  observeEvent(input$run_btn, {
    strata_res(NULL); run_err(NULL)
    run_snapshot <- design_snapshot()
    constraint_now <- input$constraint %||% "n_total"
    method_now <- input$method %||% "cobyla"
    n_starts_now <- input$n_starts %||% 8
    max_iter_now <- input$max_iter %||% 80

    # The solvers draw random initial boundaries, so an unseeded run is not
    # reproducible -- and the R Code tab promises a script that reproduces
    # what is on screen. Seed here, record the seed in the settings, and emit
    # the same set.seed() call into the generated script.
    seed_now <- suppressWarnings(as.integer(input$run_seed %||% 1L))
    if (is.na(seed_now)) seed_now <- 1L
    set.seed(seed_now)

    # ── Density (data-free) run path ─────────────────────────────────────────
    if (identical(input$data_source, "density")) {
      withProgress(message = "▶ Density optimising…", value = 0, {
        tryCatch({
          if (!"stratify_multi_density" %in% getNamespaceExports("multistrata"))
            stop("The installed multistrata predates the density mode. ",
                 "Restart R and re-source run_app.R so the updated package loads.")
          a <- gather_density_args()
          incProgress(0.1, detail = sprintf("%d strata (data-free)", prod(a$Hs)))
          t0 <- proc.time()[3]
          res <- stratify_multi_density(
            distributions = a$dists,
            H_vec         = a$Hs,
            n_total       = as.integer(input$n_total),
            N             = a$N,
            ranges        = a$ranges,
            omega         = a$om,
            constraint    = constraint_now,
            cost_total    = if (constraint_now == "cost") input$cost_total else NULL,
            cost_fixed    = if (constraint_now == "cost") input$cost_fixed else 0,
            cost_per_unit = if (constraint_now == "cost") input$cost_unit else 1,
            n_min         = as.integer(input$n_min),
            n_starts      = as.integer(n_starts_now),
            max_iter      = as.integer(max_iter_now))
          incProgress(0.9)
          strata_res(list(res = res, a = a, secs = proc.time()[3] - t0,
                          density = TRUE,
                          snapshot = run_snapshot,
                          settings = list(
                            data_source = "density",
                            n_total = input$n_total, n_min = input$n_min,
                            constraint = constraint_now,
                            cost_total = input$cost_total, cost_fixed = input$cost_fixed,
                            cost_unit = input$cost_unit,
                            method = "cobyla", wmode = "density", N = a$N,
                            n_starts = n_starts_now, max_iter = max_iter_now,
                            seed = seed_now)))
          bslib::nav_select("main_tabs", "design", session)
        }, error = function(e) run_err(conditionMessage(e)))
      })
      return()
    }

    withProgress(message = "▶ MSAO optimising…", value = 0, {
      tryCatch({
        a <- gather_args()
        incProgress(0.1, detail = sprintf("%d strata · %d starts",
                                          prod(a$Hs), n_starts_now))
        t0 <- proc.time()[3]
        res <- stratify_multi(
          data          = a$dat,
          study_vars    = a$sv,
          H_vec         = a$Hs,
          n_total       = as.integer(input$n_total),
          omega         = a$om,
          constraint    = constraint_now,
          cost_total    = if (constraint_now == "cost") input$cost_total else NULL,
          cost_fixed    = if (constraint_now == "cost") input$cost_fixed else 0,
          cost_per_unit = if (constraint_now == "cost") input$cost_unit else 1,
          cv_targets    = a$cvt,
          n_min         = as.integer(input$n_min),
          n_starts      = as.integer(n_starts_now),
          max_iter      = as.integer(max_iter_now),
          method        = method_now,
          distributions = a$dists)
        incProgress(0.9)
        strata_res(list(res = res, a = a, secs = proc.time()[3] - t0,
                        snapshot = run_snapshot,
                        settings = list(
                          data_source = input$data_source,
                          builtin_ds = input$builtin_ds,
                          file_name = input$file_up$name %||% NULL,
                          sim_N = input$sim_N, sim_p = input$sim_p,
                          sim_q = input$sim_q, sim_scenario = input$sim_scenario,
                          sim_seed = input$sim_seed,
                          n_total = input$n_total, n_min = input$n_min,
                          constraint = constraint_now,
                          cost_total = input$cost_total, cost_fixed = input$cost_fixed,
                          cost_unit = input$cost_unit,
                          cv_targets = a$cvt,
                          method = method_now, wmode = input$wmode,
                          n_starts = n_starts_now, max_iter = max_iter_now,
                          seed = seed_now)))
        bslib::nav_select("main_tabs", "design", session)
      }, error = function(e) run_err(conditionMessage(e)))
    })
  })

  output$ui_run_status <- renderUI({
    s <- strata_res()
    if (is.null(s)) return(NULL)
    div(class = "ok-box",
        sprintf("✓ Optimised %d strata in %.1f s · Φ = %s",
                nrow(s$res$stratum_stats), s$secs, fmt_num(s$res$objective, 6)))
  })

  output$ui_error <- renderUI({
    if (is.null(run_err())) return(NULL)
    div(class = "err-box", "✖ Optimisation failed:\n", run_err())
  })

  output$ui_stale <- renderUI({
    if (!isTRUE(results_stale())) return(NULL)
    div(class = "warn-box",
        tags$strong("Settings changed since this design was run. "),
        "The displayed results remain the previous completed design. Press Run Joint ",
        "Optimisation to refresh results before using follow-on analyses or downloads.")
  })

  # ── Metric strip ───────────────────────────────────────────────────────────
  output$ui_metrics <- renderUI({
    s <- strata_res(); if (is.null(s)) return(NULL)
    r <- s$res
    div(class = "metrics-row",
      div(class = "metric-card",
        div(class = "metric-value", fmt_num(r$objective, 6)),
        div(class = "metric-label", "Composite Φ")),
      div(class = "metric-card",
        div(class = "metric-value", sprintf("%.2f%%", 100 * r$optimality_gap)),
        div(class = "metric-label", "Optimality gap")),
      div(class = "metric-card",
        div(class = "metric-value",
            if (is.null(r$kkt_residual)) "n/a"
            else formatC(r$kkt_residual, format = "e", digits = 1)),
        div(class = "metric-label", "KKT residual")),
      div(class = "metric-card",
        div(class = "metric-value", nrow(r$stratum_stats)),
        div(class = "metric-label", "Strata")),
      div(class = "metric-card",
        div(class = "metric-value", sum(r$allocation)),
        div(class = "metric-label", "Sample n")),
      div(class = "metric-card",
        div(class = "metric-value",
            paste0(toupper(substr(r$.meta$mode, 1, 5)), " · ",
                   toupper(r$.meta$method))),
        div(class = "metric-label", "Mode · solver"))
    )
  })

  # ============================================================================
  # TAB 1: Introduction
  # ============================================================================
  output$ui_intro <- renderUI({
    tagList(
      div(class = "intro-section",
        div(class = "intro-section-title", div(class = "num-badge", "1"),
            "The Joint Multivariate Design Problem"),
        div(class = "intro-body",
          tags$p("Classical survey design fixes strata first and allocates the",
            " sample afterwards. multistrata solves both problems",
            tags$strong(" simultaneously"), " across several stratification",
            " variables (X) and several study variables (Y), minimising the",
            " weighted composite variance:"),
          div(class = "intro-formula",
            "Φ(B, n) = Σₖ ωₖ Vₖ(B, n), ",
            "Vₖ = Σₕ Wₕ² σ²ₕₖ (1/nₕ − 1/Nₕ)"),
          tags$p("The joint optimiser is designed to avoid the common",
            " sequential shortcut of fixing boundaries first and allocating",
            " afterwards. The Benchmark tab compares the joint design with",
            " simpler alternatives on your data.")
        )),
      div(class = "intro-section",
        div(class = "intro-section-title", div(class = "num-badge", "2"),
            "Two Weight Modes"),
        div(class = "intro-body",
          tags$p(tags$strong("Empirical:"), " stratum weights Wₕ are sample",
            " proportions, with no assumptions, but the objective is a step function."),
          tags$p(tags$strong("Parametric:"), " each Xⱼ is assigned one of 13",
            " distribution families and Wₕ is computed analytically from",
            " products of marginal CDF increments. This is the multivariate",
            " version of the density-based stratifyR approach. The objective becomes",
            " smooth, enabling the gradient-based SLSQP solver via the",
            " closed-form weight gradient ∂Wₕ/∂b = ±fⱼ(b)·ΠΔF.")
        )),
      div(class = "intro-section",
        div(class = "intro-section-title", div(class = "num-badge", "3"),
            "The MSAO Algorithm"),
        div(class = "intro-body",
          tags$p("Multi-Start Alternating Optimisation: from each random start,",
            " a boundary update (COBYLA or SLSQP) alternates with an exact",
            " generalised Neyman allocation (KKT water-filling:",
            " nₕ ∝ Aₕ = Wₕ√(Σₖωₖσ²ₕₖ)),",
            " with monotone descent. Constraints: fixed n (C1), budget",
            " c₀ + Σcₕnₕ ≤ C* (C2), CV targets per study",
            " variable (C3, convex sub-problem), minimum stratum sizes (C4)."),
          tags$p("Every returned design carries two quality diagnostics:",
            " the Cauchy-Schwarz ", tags$strong("optimality gap"),
            ", which concerns the allocation, and a scaled ",
            tags$strong("KKT residual"), " reported in two parts, one for the",
            " allocation and one for the boundaries, so allocation optimality",
            " is kept separate from any unresolved boundary-search error.")
        )),
      div(class = "intro-section",
        div(class = "intro-section-title", div(class = "num-badge", "4"),
            "How to Use This App"),
        div(class = "intro-body",
          tags$ol(
            tags$li(tags$strong("Data"), ": pick a built-in population, upload",
              " CSV/Excel, or simulate one; assign X and Y roles. Or choose",
              tags$strong(" Densities only"), " to design directly from the",
              " marginal distributions with no data."),
            tags$li(tags$strong("Design"), ": strata per variable Hⱼ, sample",
              " size n, priority weights ω."),
            tags$li(tags$strong("Weights"), ": empirical, or parametric with",
              " per-variable families (one-click MLE prefill)."),
            tags$li(tags$strong("Constraint & solver"), ": C1/C2/C3, COBYLA or",
              " SLSQP, number of multi-starts."),
            tags$li(tags$strong("Run"), ": start the optimiser, then explore the interactive boundary",
              " maps, 3-D views, density contours, scatter matrices, benchmark",
              " comparisons, CV dashboards, the ω trade-off explorer, and a",
              " reproducible R script.")))),
      div(class = "intro-section",
        div(class = "intro-section-title", div(class = "num-badge", "?"),
            "Key Terms in Plain Language"),
        div(class = "intro-body",
          tags$ul(
            tags$li(tags$strong("Strata: "), "groups of similar population units. For example, small farms, medium farms and large farms."),
            tags$li(tags$strong("Allocation: "), "how many survey units to sample from each stratum."),
            tags$li(tags$strong("CV: "), "coefficient of variation, a percentage measure of survey precision. Smaller is usually better."),
            tags$li(tags$strong("ω / omega: "), "priority weights. A larger weight tells the optimiser to care more about that study variable."),
            tags$li(tags$strong("KKT residual: "), "a technical quality check for the optimisation. Smaller values mean the solver is closer to a stable optimum."),
            tags$li(tags$strong("Optimality gap: "), "another quality check. A small gap means the design is close to the best-known lower bound.")))),
      div(class = "intro-ref",
        div(class = "intro-ref-title", "References and background"),
        tags$p("This companion app accompanies ongoing multistrata methodology work. Related foundations include:"),
        tags$ul(
          tags$li("Reddy, K.G. & Khan, M.G.M. (2020). stratifyR: An R package for optimal stratification and sample allocation for univariate populations. ", tags$em("ANZJS"), " 62(3), 383-405."),
          tags$li("Bethel, J. (1989). Sample allocation in multivariate surveys. ", tags$em("Survey Methodology"), " 15, 47-57."),
          tags$li("Neyman, J. (1934). ", tags$em("JRSS"), " 97, 558-625."),
          tags$li("Cochran, W.G. (1977). ", tags$em("Sampling Techniques"), ", 3rd ed. Wiley.")))
    )
  })

  # ============================================================================
  # TAB 2: Optimal Design
  # ============================================================================
  output$ui_design <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run an optimisation first. This tab will then show the strata boundaries, sample allocation, and downloadable tables.",
      title = "Run first to see the design",
      icon = "\U0001F9ED"))
    r <- s$res
    tagList(
      local({
        alloc <- as.integer(r$allocation %||% integer())
        top <- if (length(alloc)) which.max(alloc) else NA_integer_
        top_lab <- if (is.finite(top) && "stratum" %in% names(r$stratum_stats)) {
          as.character(r$stratum_stats$stratum[top])
        } else if (is.finite(top)) {
          paste0("stratum ", top)
        } else {
          "the largest stratum"
        }
        div(class = "result-insight",
          sprintf("Quick read: %d strata, %s sample units. The largest allocation is %s (n = %s).",
                  nrow(r$stratum_stats), format(sum(alloc), big.mark = ","),
                  top_lab, if (length(alloc)) format(max(alloc), big.mark = ",") else "n/a"))
      }),
      div(class = "explain-card",
        div(class = "explain-title", "What you are looking at"),
        tags$p("The boundaries define the groups, or strata. The allocation table shows how many units to sample from each group. Larger allocations usually go to strata that are bigger, more variable, or more important for your chosen study variables.")),
      tags$h5("Optimal strata boundaries"),
      lapply(seq_along(s$a$Xv), function(j) {
        b <- r$boundaries[[j]]
        div(class = "bounds-card",
          tags$span(class = "bounds-var", s$a$Xv[j]), ":  [",
          fmt_num(r$.meta$alpha[j], 3), " | ",
          paste(fmt_num(b, 4), collapse = " | "), " | ",
          fmt_num(r$.meta$beta[j], 3), "]")
      }),
      if (!is.null(r$cv_achieved))
        div(class = if (all(r$cv_achieved <= (r$.meta$cv_targets %||% Inf) + 1e-6))
                      "ok-box" else "warn-box",
          tags$strong("CV constraint (C3): "),
          paste(sprintf("%s = %.2f%% (target %.2f%%)", s$a$Yv,
                        100 * r$cv_achieved, 100 * r$.meta$cv_targets),
                collapse = " · ")),
      if (isTRUE(s$a$auto)) tagList(
        div(class = "ok-box",
          tags$strong("Auto-parametric selection: "),
          paste(sprintf("%s ~ %s", s$a$Xv,
                        vapply(s$a$dists, function(d) d$family, "")),
                collapse = "  ·  ")),
        tags$details(
          tags$summary(style = "cursor:pointer; font-size:0.82rem; color:#24527A;",
                       "Show AIC fit table"),
          DTOutput("dt_autofit"))),
      tags$h5("Stratum-level statistics"),
      DTOutput("dt_strata"),
      div(class = "dl-row",
        downloadButton("dl_strata", "Stratum table (CSV)", class = "btn-sm"),
        downloadButton("dl_alloc", "Allocation (CSV)", class = "btn-sm"),
        downloadButton("dl_bounds", "Boundaries (CSV)", class = "btn-sm"))
    )
  })

  output$dt_strata <- renderDT({
    s  <- strata_res(); req(s)
    # One table with everything: per-variable boundary ranges + N_h (population),
    # W_h (weight), n_h (sample), study-variable moments and A_h.
    df <- .stratum_table(s)
    round_cols <- setdiff(names(df)[vapply(df, is.numeric, TRUE)],
                          c("N_h (population)", "n_h (sample)"))
    datatable(df, rownames = FALSE,
              options = list(pageLength = 12, scrollX = TRUE, dom = "tip")) |>
      formatRound(round_cols, 4)
  })

  output$dt_autofit <- renderDT({
    s <- strata_res(); req(s, s$a$auto_tbl)
    datatable(s$a$auto_tbl, rownames = FALSE,
              options = list(pageLength = 12, dom = "tp")) |>
      formatRound(c("logLik", "AIC"), 1)
  })

  output$dl_strata <- downloadHandler(
    filename = function() "multistrata_stratum_stats.csv",
    content  = function(f) {
      s <- strata_res(); req(s)
      utils::write.csv(.stratum_table(s), f, row.names = FALSE)
    })

  output$dl_alloc <- downloadHandler(
    filename = function() "multistrata_allocation.csv",
    content  = function(f) {
      s <- strata_res(); req(s)
      utils::write.csv(.allocation_table(s), f, row.names = FALSE)
    })

  output$dl_bounds <- downloadHandler(
    filename = function() "multistrata_boundaries.csv",
    content  = function(f) {
      s <- strata_res(); req(s)
      utils::write.csv(.boundaries_table(s), f, row.names = FALSE)
    })

  # ============================================================================
  # TAB 3: Visualise (interactive, multivariable)
  # ============================================================================
  output$ui_visualise <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run an optimisation first, then this tab will show boundary maps, allocation plots, and interactive figures.",
      title = "Run first to visualise",
      icon = "\U0001F5FA"))
    if (isTRUE(s$density)) {
      p <- length(s$a$Xv)
      marg_panel <- nav_panel("Marginal densities",
        div(class = "explain-card",
          div(class = "explain-title", "Why this matters"),
          tags$p("Visual checks help you see whether the strata make practical sense. A good design should split the population into interpretable groups and place more sample effort where it improves precision.")),
        div(class = "mini-note",
            "Each stratification variable's marginal density with the optimal",
            " cut-points; shaded bands are the strata. The design is computed",
            " analytically from these curves, there is no data cloud."),
        plotOutput("pl_density", height = "640px"))
      alloc_panel <- nav_panel("Allocation",
        fluidRow(
          column(4, selectInput("alloc_metric", "Metric",
            c("Sample size nₕ" = "n_h", "Weight Wₕ" = "W_h",
              "Aₕ = Wₕ√(Σωσ²)" = "A_h"))),
          column(4, checkboxInput("alloc_sort", "Sort descending", FALSE))),
        if (has_plotly) plotly_output("pl_alloc", height = "420px")
        else plotOutput("plot_pkg_alloc", height = "420px"),
        if (has_plotly) plotly_output("pl_bubble", height = "420px"))
      wts_panel <- nav_panel("Weights", plotOutput("plot_pkg_wts", height = "420px"))

      if (p < 2 || !has_plotly)
        return(div(class = "visualise-panel", navset_pill(marg_panel, alloc_panel, wts_panel)))

      return(div(class = "visualise-panel", navset_pill(
        nav_panel("Boundary Map",
          div(class = "mini-note",
              "The joint density f(X₁,X₂)=∏ⱼ fⱼ (independent margins) with the",
              " optimal boundary grid. Cells are the strata; labels show nₕ, Wₕ."),
          div(class = "vis-controls",
            if (p > 2) uiOutput("ui_den_pair_map"),
            checkboxInput("den_map_labels", "Cell labels (nₕ, Wₕ)", value = (p == 2)),
            sliderInput("den_grid", "Grid resolution", 40, 140, 80,
                        step = 10, ticks = FALSE, width = "240px")),
          div(class = "vis-plot", plotly_output("pl_den_map", height = "620px"))),
        nav_panel("3-D Density",
          div(class = "vis-controls",
            if (p > 2) uiOutput("ui_den_pair_3d"),
            checkboxInput("den_planes", "Boundary walls", TRUE)),
          div(class = "vis-plot", plotly_output("pl_den_3d", height = "640px"))),
        nav_panel("Density Contour",
          div(class = "vis-controls", if (p > 2) uiOutput("ui_den_pair_ct")),
          div(class = "vis-plot", plotly_output("pl_den_contour", height = "600px"))),
        marg_panel, alloc_panel, wts_panel)))
    }
    if (!has_plotly)
      return(tagList(
        div(class = "warn-box",
            "Package 'plotly' not installed, showing the static package plots",
            " only. Run install.packages('plotly') for the interactive views."),
        uiOutput("ui_vis_static")))
    p <- length(s$a$Xv)
    if (p < 2)
      return(div(class = "visualise-panel", navset_pill(
        nav_panel("Allocation",
          div(class = "vis-controls",
            selectInput("alloc_metric", "Metric",
              c("Sample size n" = "n_h", "Weight W" = "W_h", "A_h" = "A_h"),
              width = "240px"),
            checkboxInput("alloc_sort", "Sort descending", FALSE)),
          div(class = "vis-plot", plotly_output("pl_alloc", height = "420px")),
          div(class = "vis-plot", plotly_output("pl_bubble", height = "420px"))),
        nav_panel("Static (package)", uiOutput("ui_vis_static")))))
    div(class = "visualise-panel", navset_pill(
      nav_panel("Boundary Map",
        div(class = "vis-controls",
          if (p > 2) uiOutput("ui_pair_sel"),
          checkboxInput("map_labels", "Cell labels (nₕ, Wₕ)", value = (p == 2)),
          sliderInput("map_pts", "Points shown", 500, 8000, 3000,
                      step = 500, ticks = FALSE, width = "280px")),
        div(class = "vis-plot", plotly_output("pl_map", height = "620px"))),
      nav_panel("3-D Explorer",
        div(class = "vis-controls",
          if (p > 2) uiOutput("ui_pair_sel3d"),
          uiOutput("ui_zvar"),
          sliderInput("pts3d", "Points shown", 500, 6000, 2500,
                      step = 500, ticks = FALSE, width = "260px"),
          checkboxInput("planes3d", "Boundary planes", TRUE)),
        div(class = "vis-plot", plotly_output("pl_3d", height = "640px"))),
      nav_panel("Density Contour",
        div(class = "vis-controls",
          if (p > 2) uiOutput("ui_pair_selkd"),
          checkboxInput("kd_pts", "Overlay points", FALSE)),
        div(class = "vis-plot", plotly_output("pl_kde", height = "600px"))),
      nav_panel("3-D Density Ridge",
        div(class = "vis-controls",
          radioButtons("ridge_view", NULL,
            c("X variables + boundaries" = "x",
              "Y variable, by stratum"   = "y"), selected = "x"),
          conditionalPanel("input.ridge_view == 'y'", uiOutput("ui_ridge_y")),
          conditionalPanel("input.ridge_view == 'x'",
            if (has_fitverse)
              actionButton("fit_fv", "⚙ Fit best density (FitVerse)",
                           class = "btn-mle", style = "margin-top:0.3rem;")
            else div(class = "mini-note",
                     "Install FitVerse for best-fit density overlays.")),
          uiOutput("ui_ridge_status")),
        div(class = "mini-note",
            "Each variable (or stratum) occupies one lane; curves are",
            " kernel densities scaled to a common height, values rescaled",
            " to [0, 1] per lane so different measurement scales share one",
            " axis. Gold walls mark the optimal boundary cuts. Drag to",
            " rotate; hover for real (unscaled) values."),
        div(class = "vis-plot", plotly_output("pl_ridge", height = "660px"))),
      nav_panel("Scatter Matrix",
        div(class = "mini-note",
            "All stratification variables (and optionally study variables),",
            " coloured by optimal stratum. Drag to select, double-click to reset."),
        div(class = "vis-controls", checkboxInput("splom_y", "Include study variables", FALSE)),
        div(class = "vis-plot", plotly_output("pl_splom", height = "680px"))),
      nav_panel("Allocation",
        div(class = "vis-controls",
          selectInput("alloc_metric", "Metric",
            c("Sample size nₕ" = "n_h", "Weight Wₕ" = "W_h",
              "Aₕ = Wₕ√(Σωσ²)" = "A_h"), width = "260px"),
          checkboxInput("alloc_sort", "Sort descending", FALSE)),
        div(class = "vis-plot", plotly_output("pl_alloc", height = "420px")),
        tags$h6(style = "margin-top:1rem;", "Allocation structure: nₕ vs Wₕ",
                tags$span(class = "mini-note",
                          "  (bubble size = Aₕ; Neyman pulls samples toward",
                          " high-Aₕ strata, away from proportionality)")),
        div(class = "vis-plot", plotly_output("pl_bubble", height = "420px"))),
      nav_panel("Static (package)", uiOutput("ui_vis_static"))
    ))
  })

  # pair selectors for p > 2
  pair_ui <- function(id) {
    s <- strata_res(); req(s, !isTRUE(s$density))
    selectizeInput(id, "Variable pair", choices = s$a$Xv,
                   selected = s$a$Xv[1:2], multiple = TRUE,
                   options = list(maxItems = 2))
  }
  output$ui_pair_sel   <- renderUI(pair_ui("map_pair"))
  output$ui_pair_sel3d <- renderUI(pair_ui("pair3d"))
  output$ui_pair_selkd <- renderUI(pair_ui("pairkd"))

  output$ui_zvar <- renderUI({
    s <- strata_res(); req(s, !isTRUE(s$density))
    selectInput("zvar", "Z axis (study variable)", choices = s$a$Yv)
  })

  get_pair <- function(sel) {
    s <- strata_res(); req(s)
    Xv <- s$a$Xv
    if (length(Xv) <= 2) return(c(1, 2)[seq_len(min(2, length(Xv)))])
    if (is.null(sel) || length(sel) < 2) return(c(1, 2))
    match(sel[1:2], Xv)
  }

  vis_data <- reactive({
    s <- strata_res(); req(s, !isTRUE(s$density))
    sid <- attr(s$res$stratum_stats, "strat_id")
    data.frame(s$a$dat, s$a$sv, .stratum = factor(sid),
               check.names = FALSE)
  })

  # ── Interactive 2-D boundary map ───────────────────────────────────────────
  output$pl_map <- render_plotly({
    s <- strata_res(); req(s, !isTRUE(s$density), has_plotly)
    jj <- get_pair(input$map_pair); req(all(!is.na(jj)))
    v1 <- s$a$Xv[jj[1]]; v2 <- s$a$Xv[jj[2]]
    df <- vis_data()
    n_show <- min(nrow(df), input$map_pts %||% 3000)
    set.seed(1); idx <- sample.int(nrow(df), n_show)
    d <- df[idx, ]
    r  <- s$res
    b1 <- r$boundaries[[jj[1]]]; b2 <- r$boundaries[[jj[2]]]
    rx <- c(r$.meta$alpha[jj[1]], r$.meta$beta[jj[1]])
    ry <- c(r$.meta$alpha[jj[2]], r$.meta$beta[jj[2]])

    shp <- c(
      lapply(b1, function(b) list(type = "line", x0 = b, x1 = b,
        y0 = ry[1], y1 = ry[2],
        line = list(color = "#0D1B2A", width = 2, dash = "solid"))),
      lapply(b2, function(b) list(type = "line", x0 = rx[1], x1 = rx[2],
        y0 = b, y1 = b,
        line = list(color = "#0D1B2A", width = 2, dash = "solid"))))

    ann <- list()
    if (isTRUE(input$map_labels) && length(s$a$Xv) == 2) {
      ss <- r$stratum_stats
      H1 <- r$.meta$H_vec[1]
      cells <- cell_rects(b1, b2, rx, ry)
      ann <- lapply(cells, function(cl) {
        h <- (cl$iy - 1L) * H1 + cl$ix        # mixed-radix: first var fastest
        list(x = cl$cx, y = cl$cy,
             text = sprintf("<b>%s</b><br>n=%d · W=%.3f",
                            ss$stratum[h], ss$n_h[h], ss$W_h[h]),
             showarrow = FALSE,
             font = list(size = 10, color = "#0D1B2A"),
             bgcolor = "rgba(255,255,255,0.75)", borderpad = 2)
      })
    }

    plot_ly(d, x = ~get(v1), y = ~get(v2), color = ~.stratum,
            colors = pal_strata[seq_len(nlevels(d$.stratum))],
            type = "scattergl", mode = "markers",
            marker = list(size = 5, opacity = 0.55),
            text = ~paste0("stratum ", .stratum),
            hovertemplate = paste0(v1, ": %{x:.3f}<br>", v2,
                                   ": %{y:.3f}<br>%{text}<extra></extra>")) |>
      layout(shapes = shp, annotations = ann,
             xaxis = list(title = v1, zeroline = FALSE),
             yaxis = list(title = v2, zeroline = FALSE),
             legend = list(title = list(text = "Stratum"), orientation = "h",
                           y = -0.15),
             margin = list(t = 20))
  })

  # ── 3-D explorer with boundary planes ──────────────────────────────────────
  output$pl_3d <- render_plotly({
    s <- strata_res(); req(s, !isTRUE(s$density), has_plotly, input$zvar)
    jj <- get_pair(input$pair3d); req(all(!is.na(jj)))
    v1 <- s$a$Xv[jj[1]]; v2 <- s$a$Xv[jj[2]]; vz <- input$zvar
    df <- vis_data()
    n_show <- min(nrow(df), input$pts3d %||% 2500)
    set.seed(2); d <- df[sample.int(nrow(df), n_show), ]
    r  <- s$res
    b1 <- r$boundaries[[jj[1]]]; b2 <- r$boundaries[[jj[2]]]
    rx <- range(d[[v1]]); ry <- range(d[[v2]]); rz <- range(d[[vz]])

    fig <- plot_ly() |>
      add_trace(data = d, x = ~get(v1), y = ~get(v2), z = ~get(vz),
                color = ~.stratum,
                colors = pal_strata[seq_len(nlevels(d$.stratum))],
                type = "scatter3d", mode = "markers",
                marker = list(size = 2.6, opacity = 0.65),
                text = ~paste0("stratum ", .stratum),
                hovertemplate = paste0(v1, ": %{x:.3f}<br>", v2, ": %{y:.3f}<br>",
                                       vz, ": %{z:.3f}<br>%{text}<extra></extra>"))
    if (isTRUE(input$planes3d)) {
      teal <- "rgba(42,157,143,0.9)"
      for (b in b1)
        fig <- fig |> add_surface(
          x = matrix(b, 2, 2), y = matrix(c(ry[1], ry[2], ry[1], ry[2]), 2, 2),
          z = matrix(c(rz[1], rz[1], rz[2], rz[2]), 2, 2),
          opacity = 0.22, showscale = FALSE,
          colorscale = list(c(0, teal), c(1, teal)),
          hoverinfo = "skip")
      for (b in b2)
        fig <- fig |> add_surface(
          x = matrix(c(rx[1], rx[2], rx[1], rx[2]), 2, 2), y = matrix(b, 2, 2),
          z = matrix(c(rz[1], rz[1], rz[2], rz[2]), 2, 2),
          opacity = 0.22, showscale = FALSE,
          colorscale = list(c(0, teal), c(1, teal)),
          hoverinfo = "skip")
    }
    fig |> layout(
      scene = list(
        xaxis = list(title = v1), yaxis = list(title = v2),
        zaxis = list(title = vz),
        camera = list(eye = list(x = 1.6, y = -1.6, z = 0.9))),
      legend = list(title = list(text = "Stratum")),
      margin = list(t = 10))
  })

  # ── Density contour + boundary grid ────────────────────────────────────────
  output$pl_kde <- render_plotly({
    s <- strata_res(); req(s, !isTRUE(s$density), has_plotly, has_mass)
    jj <- get_pair(input$pairkd); req(all(!is.na(jj)))
    v1 <- s$a$Xv[jj[1]]; v2 <- s$a$Xv[jj[2]]
    df <- vis_data()
    r  <- s$res
    kd <- MASS::kde2d(df[[v1]], df[[v2]], n = 90)
    b1 <- r$boundaries[[jj[1]]]; b2 <- r$boundaries[[jj[2]]]
    ry <- range(kd$y); rx <- range(kd$x)

    shp <- c(
      lapply(b1, function(b) list(type = "line", x0 = b, x1 = b,
        y0 = ry[1], y1 = ry[2], line = list(color = "#E9A826", width = 2.5))),
      lapply(b2, function(b) list(type = "line", x0 = rx[1], x1 = rx[2],
        y0 = b, y1 = b, line = list(color = "#E9A826", width = 2.5))))

    fig <- plot_ly(x = kd$x, y = kd$y, z = t(kd$z), type = "contour",
                   colorscale = "Blues", reversescale = TRUE,
                   contours = list(showlines = FALSE), ncontours = 22,
                   colorbar = list(title = "density")) |>
      layout(shapes = shp,
             xaxis = list(title = v1), yaxis = list(title = v2),
             margin = list(t = 20))
    if (isTRUE(input$kd_pts)) {
      set.seed(3); d <- df[sample.int(nrow(df), min(nrow(df), 1500)), ]
      fig <- fig |> add_trace(data = d, x = ~get(v1), y = ~get(v2),
        type = "scattergl", mode = "markers", inherit = FALSE,
        marker = list(size = 3, color = "rgba(13,27,42,0.35)"),
        hoverinfo = "skip", showlegend = FALSE)
    }
    fig
  })

  # ── Scatter-plot matrix (any p) ────────────────────────────────────────────
  output$pl_splom <- render_plotly({
    s <- strata_res(); req(s, !isTRUE(s$density), has_plotly)
    df <- vis_data()
    vars <- s$a$Xv
    if (isTRUE(input$splom_y)) vars <- c(vars, s$a$Yv)
    set.seed(4); d <- df[sample.int(nrow(df), min(nrow(df), 2500)), ]
    dims <- lapply(vars, function(v) list(label = v, values = d[[v]]))
    plot_ly(data = d, type = "splom", dimensions = dims,
            color = ~.stratum,
            colors = pal_strata[seq_len(nlevels(d$.stratum))],
            marker = list(size = 3.5, opacity = 0.55)) |>
      layout(showlegend = TRUE,
             legend = list(title = list(text = "Stratum")),
             plot_bgcolor = "#F8FAFD", margin = list(t = 20))
  })

  # ── Allocation bars + bubble ───────────────────────────────────────────────
  output$pl_alloc <- render_plotly({
    s <- strata_res(); req(s, has_plotly)
    ss <- s$res$stratum_stats
    m  <- input$alloc_metric %||% "n_h"
    d  <- data.frame(stratum = ss$stratum, val = ss[[m]],
                     n = ss$n_h, W = ss$W_h, A = ss$A_h)
    if (isTRUE(input$alloc_sort)) d <- d[order(-d$val), ]
    d$stratum <- factor(d$stratum, levels = d$stratum)
    lab <- c(n_h = "nₕ", W_h = "Wₕ", A_h = "Aₕ")[m]
    plot_ly(d, x = ~stratum, y = ~val, type = "bar",
            marker = list(color = pal_strata[seq_len(nrow(d))]),
            text = ~sprintf("n = %d<br>W = %.3f<br>A = %.4f", n, W, A),
            textposition = "none",
            hovertemplate = paste0("<b>%{x}</b><br>", lab,
                                   " = %{y:.4f}<br>%{text}<extra></extra>")) |>
      layout(xaxis = list(title = "", tickangle = -40),
             yaxis = list(title = lab),
             margin = list(t = 15, b = 90))
  })

  output$pl_bubble <- render_plotly({
    s <- strata_res(); req(s, has_plotly)
    ss <- s$res$stratum_stats
    n_tot <- sum(ss$n_h)
    plot_ly(ss, x = ~W_h, y = ~n_h, type = "scatter", mode = "markers+text",
            name = "strata", text = ~stratum, textposition = "top center",
            textfont = list(size = 9, color = "#334155"),
            marker = list(size = ~30 * A_h / max(A_h) + 8,
                          color = pal_strata[seq_len(nrow(ss))],
                          opacity = 0.8,
                          line = list(color = "#0D1B2A", width = 1)),
            hovertemplate = paste0("%{text}<br>W=%{x:.3f} · n=%{y}",
                                   "<extra></extra>")) |>
      add_trace(x = ss$W_h, y = n_tot * ss$W_h, type = "scatter",
                mode = "lines", inherit = FALSE, showlegend = TRUE,
                name = "proportional line n = n·Wₕ",
                line = list(color = "#94A3B8", dash = "dash"),
                hoverinfo = "skip") |>
      layout(xaxis = list(title = "stratum weight Wₕ"),
             yaxis = list(title = "allocated nₕ"),
             legend = list(orientation = "h", y = -0.18),
             margin = list(t = 15))
  })

  # ── Static package plots ───────────────────────────────────────────────────
  output$ui_vis_static <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run an optimisation first, then this tab will compare the joint design with simpler benchmark designs.",
      title = "Run first to compare designs",
      icon = "\U00002696"))
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "Why this matters"),
        tags$p("Benchmarking answers a simple question: is the joint design better than easier shortcuts? Values above 1 mean the shortcut has larger variance than the joint design, so the joint design is more efficient.")),
      div(class = "mini-note",
          "Publication-quality figures produced by the package itself:",
          " plot_strata() and plot.multistrata()."),
      fluidRow(
        column(12, plotOutput("plot_pkg_strata", height = "520px"))),
      fluidRow(
        column(6, plotOutput("plot_pkg_alloc", height = "360px")),
        column(6, plotOutput("plot_pkg_wts",  height = "360px")))
    )
  })

  output$plot_pkg_strata <- renderPlot({
    s <- strata_res(); req(s)
    tryCatch(
      plot_strata(s$res, data = s$a$dat, study_vars = s$a$sv),
      error = function(e) {
        plot.new()
        title(main = paste("plot_strata() unavailable:", conditionMessage(e)),
              cex.main = 0.9)
      })
  })
  # Density-mode marginal densities with optimal cut-points
  output$pl_density <- renderPlot({
    s <- strata_res(); req(s, isTRUE(s$density))
    .plot_density_margins(s$res)
  })

  # ── Density-mode interactive joint-density views (mirror the data mode) ─────
  den_pair_ui <- function(id) {
    s <- strata_res(); req(s)
    selectizeInput(id, "Variable pair", choices = s$a$Xv,
                   selected = s$a$Xv[1:2], multiple = TRUE,
                   options = list(maxItems = 2))
  }
  output$ui_den_pair_map <- renderUI(den_pair_ui("den_pair_map"))
  output$ui_den_pair_3d  <- renderUI(den_pair_ui("den_pair_3d"))
  output$ui_den_pair_ct  <- renderUI(den_pair_ui("den_pair_ct"))

  den_get_pair <- function(sel) {
    s <- strata_res(); req(s); Xv <- s$a$Xv
    if (length(Xv) <= 2) return(seq_len(min(2L, length(Xv))))
    if (is.null(sel) || length(sel) < 2) return(c(1L, 2L))
    match(sel[1:2], Xv)
  }

  # Analytic joint density of a variable pair on a grid (margins independent).
  den_grid <- function(s, jj, n = 80L) {
    m  <- s$res$.meta
    bp <- get(".build_pdf", asNamespace("multistrata"))
    a1 <- m$alpha[jj[1]]; b1 <- m$beta[jj[1]]
    a2 <- m$alpha[jj[2]]; b2 <- m$beta[jj[2]]
    gx <- seq(a1, b1, length.out = n); gy <- seq(a2, b2, length.out = n)
    f1 <- bp(m$distributions[[jj[1]]], a1, b1)
    f2 <- bp(m$distributions[[jj[2]]], a2, b2)
    d1 <- vapply(gx, function(x) tryCatch(as.numeric(f1(x)), error = function(e) 0), 1)
    d2 <- vapply(gy, function(y) tryCatch(as.numeric(f2(y)), error = function(e) 0), 1)
    d1[!is.finite(d1)] <- 0; d2[!is.finite(d2)] <- 0
    list(gx = gx, gy = gy, z = outer(d2, d1),          # z[y, x]
         b1 = s$res$boundaries[[jj[1]]], b2 = s$res$boundaries[[jj[2]]],
         rx = c(a1, b1), ry = c(a2, b2),
         v1 = m$x_names[jj[1]], v2 = m$x_names[jj[2]])
  }

  den_bound_shapes <- function(g) {
    c(lapply(g$b1, function(b) list(type = "line", x0 = b, x1 = b,
        y0 = g$ry[1], y1 = g$ry[2], line = list(color = "#0D1B2A", width = 2))),
      lapply(g$b2, function(b) list(type = "line", x0 = g$rx[1], x1 = g$rx[2],
        y0 = b, y1 = b, line = list(color = "#0D1B2A", width = 2))))
  }

  # Boundary map: joint-density heatmap + boundary grid + cell labels
  output$pl_den_map <- render_plotly({
    s <- strata_res(); req(s, isTRUE(s$density), has_plotly)
    jj <- den_get_pair(input$den_pair_map); req(all(!is.na(jj)))
    g  <- den_grid(s, jj, as.integer(input$den_grid %||% 80))
    r  <- s$res
    ann <- list()
    if (isTRUE(input$den_map_labels) && length(s$a$Xv) == 2) {
      ss <- r$stratum_stats; H1 <- r$.meta$H_vec[1]
      ann <- lapply(cell_rects(g$b1, g$b2, g$rx, g$ry), function(cl) {
        h <- (cl$iy - 1L) * H1 + cl$ix
        list(x = cl$cx, y = cl$cy,
             text = sprintf("<b>%s</b><br>n=%d · W=%.3f",
                            ss$stratum[h], ss$n_h[h], ss$W_h[h]),
             showarrow = FALSE, font = list(size = 10, color = "#0D1B2A"),
             bgcolor = "rgba(255,255,255,0.78)", borderpad = 2)
      })
    }
    plot_ly(x = g$gx, y = g$gy, z = g$z, type = "heatmap",
            colorscale = "YlGnBu", reversescale = TRUE,
            hovertemplate = paste0(g$v1, ": %{x:.3f}<br>", g$v2,
              ": %{y:.3f}<br>density: %{z:.4f}<extra></extra>")) |>
      layout(shapes = den_bound_shapes(g), annotations = ann,
             xaxis = list(title = g$v1), yaxis = list(title = g$v2),
             title = list(text = sprintf("Joint density f(%s, %s) + optimal boundaries",
                                         g$v1, g$v2), font = list(size = 13)))
  })

  # Density contour + boundary grid
  output$pl_den_contour <- render_plotly({
    s <- strata_res(); req(s, isTRUE(s$density), has_plotly)
    jj <- den_get_pair(input$den_pair_ct); req(all(!is.na(jj)))
    g  <- den_grid(s, jj, 90L)
    plot_ly(x = g$gx, y = g$gy, z = g$z, type = "contour",
            colorscale = "YlGnBu", reversescale = TRUE,
            contours = list(coloring = "heatmap")) |>
      layout(shapes = den_bound_shapes(g),
             xaxis = list(title = g$v1), yaxis = list(title = g$v2),
             title = list(text = sprintf("Density contours: f(%s, %s)",
                                         g$v1, g$v2), font = list(size = 13)))
  })

  # 3-D joint-density surface with optional boundary walls
  output$pl_den_3d <- render_plotly({
    s <- strata_res(); req(s, isTRUE(s$density), has_plotly)
    jj <- den_get_pair(input$den_pair_3d); req(all(!is.na(jj)))
    g  <- den_grid(s, jj, 70L)
    zmax <- max(g$z, 1e-9)
    wall <- list(list(0, "#B8860B"), list(1, "#B8860B"))
    fig <- plot_ly() |>
      add_surface(x = g$gx, y = g$gy, z = g$z, colorscale = "YlGnBu",
                  reversescale = TRUE, opacity = 0.92,
                  colorbar = list(title = "f"))
    if (isTRUE(input$den_planes %||% TRUE)) {
      for (b in g$b1)
        fig <- add_surface(fig, x = matrix(b, 2, 2),
          y = matrix(c(g$ry[1], g$ry[1], g$ry[2], g$ry[2]), 2, 2),
          z = matrix(c(0, zmax, 0, zmax), 2, 2), showscale = FALSE,
          opacity = 0.25, colorscale = wall, hoverinfo = "skip")
      for (b in g$b2)
        fig <- add_surface(fig,
          x = matrix(c(g$rx[1], g$rx[1], g$rx[2], g$rx[2]), 2, 2),
          y = matrix(b, 2, 2), z = matrix(c(0, zmax, 0, zmax), 2, 2),
          showscale = FALSE, opacity = 0.25, colorscale = wall, hoverinfo = "skip")
    }
    fig |> layout(scene = list(
      xaxis = list(title = g$v1), yaxis = list(title = g$v2),
      zaxis = list(title = "density"),
      camera = list(eye = list(x = 1.8, y = -1.6, z = 0.8)),
      aspectratio = list(x = 1.2, y = 1.2, z = 0.7)))
  })

  output$plot_pkg_alloc <- renderPlot({
    s <- strata_res(); req(s); plot(s$res, which = "allocation")
  })
  output$plot_pkg_wts <- renderPlot({
    s <- strata_res(); req(s); plot(s$res, which = "weights")
  })

  # ── 3-D density ridge (X + boundaries, or Y by stratum) ───────────────────
  output$ui_ridge_y <- renderUI({
    s <- strata_res(); req(s)
    selectInput("ridge_yvar", NULL, choices = s$a$Yv)
  })

  fv_fits <- reactiveVal(NULL)

  observeEvent(input$fit_fv, {
    s <- strata_res(); req(s, has_fitverse)
    withProgress(message = "FitVerse: selecting best families…", value = 0, {
      fits <- lapply(seq_along(s$a$Xv), function(j) {
        incProgress(1 / length(s$a$Xv), detail = s$a$Xv[j])
        x <- s$a$dat[[j]]; x <- x[is.finite(x)]
        if (length(x) > 4000) { set.seed(9); x <- sample(x, 4000) }
        tryCatch({
          ft <- FitVerse::fitverse(x, plot = FALSE, verbose = FALSE)
          nm <- tryCatch(
            FitVerse::get_dist_entry(ft$best_fit$dist_key)$display_name,
            error = function(e) ft$best_fit$dist_key)
          list(fit = ft, label = nm)
        }, error = function(e) NULL)
      })
      names(fits) <- s$a$Xv
      fv_fits(fits)
    })
    showNotification("FitVerse fits ready, overlaid on the ridge.",
                     type = "message")
  })

  # reset fits when a new optimisation replaces the data
  observeEvent(strata_res(), fv_fits(NULL))

  output$ui_ridge_status <- renderUI({
    fv <- fv_fits()
    if (is.null(fv)) return(NULL)
    ok <- !vapply(fv, is.null, TRUE)
    div(class = "ok-box", style = "margin-top:0.2rem;",
        paste0("FitVerse: ",
               paste(sprintf("%s → %s", names(fv)[ok],
                             vapply(fv[ok], function(f) f$label, "")),
                     collapse = " · ")))
  })

  output$pl_ridge <- render_plotly({
    s <- strata_res(); req(s, has_plotly)
    r <- s$res
    fig <- plot_ly()

    hex2rgba <- function(hex, a) {
      v <- grDevices::col2rgb(hex)
      sprintf("rgba(%d,%d,%d,%.2f)", v[1], v[2], v[3], a)
    }

    if (identical(input$ridge_view %||% "x", "x")) {
      # ── lanes = stratification variables, gold walls = boundaries ─────────
      Xv <- s$a$Xv
      fv <- fv_fits()
      lane_names <- character(length(Xv))
      for (j in seq_along(Xv)) {
        x  <- s$a$dat[[j]]; x <- x[is.finite(x)]
        a_j <- r$.meta$alpha[j]; b_j <- r$.meta$beta[j]
        kd <- stats::density(x, n = 220, from = a_j, to = b_j)
        gx <- (kd$x - a_j) / (b_j - a_j)
        zmax <- max(kd$y)
        gz <- kd$y / zmax
        col <- pal_strata[j]
        lane_names[j] <- Xv[j]

        # density curtain (vertical ribbon under the curve)
        fig <- fig |> add_surface(
          x = rbind(gx, gx), y = matrix(j, 2, length(gx)),
          z = rbind(rep(0, length(gx)), gz),
          opacity = 0.55, showscale = FALSE,
          colorscale = list(c(0, hex2rgba(col, 0.85)),
                            c(1, hex2rgba(col, 0.85))),
          hoverinfo = "skip")
        # crest line with real-value hover
        fig <- fig |> add_trace(
          x = gx, y = rep(j, length(gx)), z = gz,
          type = "scatter3d", mode = "lines", name = Xv[j],
          line = list(color = col, width = 5),
          customdata = kd$x,
          hovertemplate = paste0(Xv[j], " = %{customdata:.3f}",
                                 "<br>density = %{z:.3f}<extra>",
                                 Xv[j], "</extra>"))
        # FitVerse best-fit overlay
        if (!is.null(fv) && !is.null(fv[[Xv[j]]])) {
          fz <- tryCatch(FitVerse::dfitverse(kd$x, fv[[Xv[j]]]$fit),
                         error = function(e) NULL)
          if (!is.null(fz) && all(is.finite(fz))) {
            fig <- fig |> add_trace(
              x = gx, y = rep(j, length(gx)), z = fz / zmax,
              type = "scatter3d", mode = "lines",
              name = paste0(Xv[j], ": ", fv[[Xv[j]]]$label),
              line = list(color = "#0D1B2A", width = 4, dash = "dash"),
              customdata = kd$x,
              hovertemplate = paste0("fitted ", fv[[Xv[j]]]$label,
                                     "<br>", Xv[j], " = %{customdata:.3f}",
                                     "<extra></extra>"))
          }
        }
        # boundary walls (gold)
        dens_at <- stats::approxfun(gx, gz, rule = 2)
        for (b in r$boundaries[[j]]) {
          nb <- (b - a_j) / (b_j - a_j)
          zt <- max(dens_at(nb), 0.05)
          fig <- fig |> add_trace(
            x = c(nb, nb), y = c(j, j), z = c(0, zt),
            type = "scatter3d", mode = "lines+markers",
            line = list(color = "#E9A826", width = 9),
            marker = list(size = 3.5, color = "#E9A826",
                          symbol = "diamond"),
            showlegend = FALSE,
            hovertemplate = paste0("boundary ", Xv[j], " = ",
                                   fmt_num(b, 4), "<extra></extra>"))
        }
      }
      yaxis_cfg <- list(title = "", tickvals = seq_along(Xv),
                        ticktext = lane_names, range = c(0.4, length(Xv) + 0.6))
    } else {
      # ── lanes = strata, curves = within-stratum density of chosen Y ───────
      req(input$ridge_yvar)
      df  <- vis_data()
      yv  <- input$ridge_yvar
      ss  <- r$stratum_stats
      H   <- nlevels(df$.stratum)
      rng <- range(df[[yv]], finite = TRUE)
      # common height scaling across strata
      kds <- lapply(seq_len(H), function(h) {
        yy <- df[[yv]][df$.stratum == levels(df$.stratum)[h]]
        yy <- yy[is.finite(yy)]
        if (length(yy) < 5) return(NULL)
        stats::density(yy, n = 180, from = rng[1], to = rng[2])
      })
      zmax <- max(vapply(kds, function(k) if (is.null(k)) 0 else max(k$y), 1))
      lane_names <- character(H)
      for (h in seq_len(H)) {
        lane_names[h] <- ss$stratum[h]
        kd <- kds[[h]]; if (is.null(kd)) next
        gx <- (kd$x - rng[1]) / diff(rng)
        gz <- kd$y / zmax
        col <- pal_strata[(h - 1) %% length(pal_strata) + 1]
        fig <- fig |> add_surface(
          x = rbind(gx, gx), y = matrix(h, 2, length(gx)),
          z = rbind(rep(0, length(gx)), gz),
          opacity = 0.5, showscale = FALSE,
          colorscale = list(c(0, hex2rgba(col, 0.85)),
                            c(1, hex2rgba(col, 0.85))),
          hoverinfo = "skip")
        fig <- fig |> add_trace(
          x = gx, y = rep(h, length(gx)), z = gz,
          type = "scatter3d", mode = "lines",
          name = sprintf("%s (n=%d)", ss$stratum[h], ss$n_h[h]),
          line = list(color = col, width = 5),
          customdata = kd$x,
          hovertemplate = paste0(ss$stratum[h], "<br>", yv,
                                 " = %{customdata:.3f}<br>density = %{z:.3f}",
                                 "<extra></extra>"))
      }
      yaxis_cfg <- list(title = "", tickvals = seq_len(H),
                        ticktext = lane_names, range = c(0.4, H + 0.6),
                        tickfont = list(size = 9))
    }

    fig |> layout(
      scene = list(
        xaxis = list(title = "normalised value (hover for real units)"),
        yaxis = yaxis_cfg,
        zaxis = list(title = "density (scaled)"),
        camera = list(eye = list(x = 1.9, y = -1.5, z = 0.75)),
        aspectratio = list(x = 1.4, y = 1.1, z = 0.6)),
      legend = list(orientation = "h", y = -0.05),
      margin = list(t = 10))
  })

  # ============================================================================
  # TAB 4: Benchmark designs
  # ============================================================================
  bench_res <- reactive({
    s <- strata_res(); req(s)
    validate(need(!results_stale(), "Settings changed; re-run the design to refresh this benchmark."))
    withProgress(message = "Benchmarking designs…", value = 0.4,
      if (isTRUE(s$density)) {
        benchmark_density_designs(result_joint = s$res)
      } else {
        benchmark_designs(
          data = s$a$dat, study_vars = s$a$sv,
          H_vec = s$res$.meta$H_vec,
          n_total = as.integer(sum(s$res$allocation)),
          omega = s$a$om, result_joint = s$res)
      })
  })

  output$ui_benchmark <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state())
    if (isTRUE(results_stale()))
      return(div(class = "warn-box",
        tags$strong("Settings changed. "),
        "Re-run the design to refresh the benchmark."))
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "Why this matters"),
        tags$p("Benchmarking answers a simple question: is the joint design better than easier shortcuts? Values above 1 mean the shortcut has larger variance than the joint design, so the joint design is more efficient.")),
      div(class = "mini-note",
        if (isTRUE(s$density))
          "Density-mode benchmark compares the analytic joint optimum against proportional, equal and density-quantile + Neyman alternatives. "
        else "",
        "Compares the joint optimum with simpler designs such as",
        " proportional, equal and sequential boundaries-then-allocation",
        " alternatives. RE > 1 means the joint design is more efficient.",
        " Cost-constrained designs are compared at",
        " their achieved total sample size."),
      if (has_plotly) plotly_output("pl_bench", height = "330px"),
      DTOutput("dt_bench"),
      div(class = "dl-row",
          downloadButton("dl_bench", "Benchmark (CSV)", class = "btn-sm"))
    )
  })
  outputOptions(output, "ui_benchmark", suspendWhenHidden = FALSE)

  output$pl_bench <- render_plotly({
    b <- bench_res(); req(b, has_plotly)
    d <- b[order(b$relative_efficiency), ]
    d$method <- factor(d$method, levels = d$method)
    cols <- ifelse(d$relative_efficiency <= 1 + 1e-9, "#2A9D8F", "#E15759")
    plot_ly(d, y = ~method, x = ~relative_efficiency, type = "bar",
            orientation = "h", marker = list(color = cols),
            hovertemplate = "%{y}<br>RE = %{x:.3f}<extra></extra>") |>
      layout(shapes = list(list(type = "line", x0 = 1, x1 = 1, y0 = -0.5,
             y1 = nrow(d) - 0.5, line = list(color = "#0D1B2A", dash = "dash"))),
             xaxis = list(title = "Relative efficiency vs joint (Φ / Φ_joint)"),
             yaxis = list(title = ""), margin = list(l = 210, t = 15))
  })

  output$dt_bench <- renderDT({
    b <- bench_res(); req(b)
    datatable(b, rownames = FALSE, options = list(dom = "t")) |>
      formatRound("composite_variance", 8) |>
      formatRound(c("relative_efficiency"), 4)
  })

  output$dl_bench <- downloadHandler(
    filename = function() "multistrata_benchmark.csv",
    content = function(f) utils::write.csv(bench_res(), f, row.names = FALSE))

  # ============================================================================
  # TAB 5: Precision (CV) dashboard
  # ============================================================================
  precision_data <- reactive({
    s <- strata_res(); req(s)
    ss <- s$res$stratum_stats
    q  <- length(s$a$Yv)
    s2 <- as.matrix(ss[, paste0("sigma2_k", seq_len(q)), drop = FALSE])
    mu <- as.matrix(ss[, paste0("mu_k",     seq_len(q)), drop = FALSE])
    cvs <- design_cv(s2, mu, ss$W_h, s$res$allocation, ss$N_h)
    fpc <- 1 / pmax(as.numeric(s$res$allocation), 1) - 1 / pmax(ss$N_h, 1)
    Vk  <- vapply(seq_len(q), function(k) sum(ss$W_h^2 * s2[, k] * fpc), 1)
    lb  <- variance_lower_bound(s2, ss$W_h, ss$N_h, s$a$om,
                                sum(s$res$allocation))
    list(cvs = as.numeric(cvs), Vk = Vk, lb = lb, q = q)
  })

  output$ui_precision <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run an optimisation first, then this tab will show the achieved CV for each study variable.",
      title = "Run first to check precision",
      icon = "\U0001F3AF"))
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "What CV means"),
        tags$p("CV is a percentage measure of uncertainty. Smaller CV means better precision. Use this tab to check whether each study variable is estimated accurately enough for your purpose.")),
      if (has_plotly) plotly_output("pl_cv", height = "330px"),
      tags$h6("Variance decomposition"),
      DTOutput("dt_prec"),
      uiOutput("ui_lb_note")
    )
  })

  output$pl_cv <- render_plotly({
    s <- strata_res(); pd <- precision_data(); req(s, pd, has_plotly)
    d <- data.frame(y = s$a$Yv, cv = 100 * pd$cvs)
    tg <- if (!is.null(s$res$.meta$cv_targets)) 100 * s$res$.meta$cv_targets else NULL
    fig <- plot_ly(d, x = ~cv, y = ~y, type = "bar", orientation = "h",
                   name = "achieved CV",
                   marker = list(color = "#24527A"),
                   hovertemplate = "%{y}: CV = %{x:.2f}%<extra></extra>")
    if (!is.null(tg))
      fig <- fig |> add_trace(x = tg, y = d$y, type = "scatter",
        mode = "markers", inherit = FALSE, name = "target γ",
        marker = list(symbol = "line-ns-open", size = 26,
                      color = "#E15759", line = list(width = 3)),
        hovertemplate = "target %{x:.2f}%<extra></extra>")
    fig |> layout(xaxis = list(title = "coefficient of variation (%)"),
                  yaxis = list(title = ""),
                  legend = list(orientation = "h", y = -0.25),
                  margin = list(l = 130, t = 15))
  })

  output$dt_prec <- renderDT({
    s <- strata_res(); pd <- precision_data(); req(s, pd)
    d <- data.frame(
      `study variable` = s$a$Yv,
      `omega` = round(s$a$om, 3),
      `V_k` = pd$Vk,
      `CV (%)` = round(100 * pd$cvs, 3),
      check.names = FALSE)
    datatable(d, rownames = FALSE, options = list(dom = "t")) |>
      formatSignif("V_k", 5)
  })

  output$ui_lb_note <- renderUI({
    s <- strata_res(); pd <- precision_data(); req(s, pd)
    div(class = "ok-box",
      sprintf(paste0("Composite Φ = %s · Cauchy-Schwarz lower bound = %s · ",
                     "gap = %.2f%% · KKT residual = %s ",
                     "(allocation %s, boundary %s)"),
              fmt_num(s$res$objective, 6), fmt_num(pd$lb, 6),
              100 * s$res$optimality_gap,
              formatC(s$res$kkt_residual, format = "e", digits = 2),
              if (is.null(s$res$kkt_alloc)) "n/a"
                else formatC(s$res$kkt_alloc, format = "e", digits = 2),
              if (is.null(s$res$kkt_boundary)) "n/a"
                else formatC(s$res$kkt_boundary, format = "e", digits = 2)))
  })

  # ============================================================================
  # TAB 6: ω Explorer (priority trade-off sweep)
  # ============================================================================
  output$ui_omega_explorer <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run an optimisation first, then explore how changing priority weights ω shifts precision between variables.",
      title = "Run first to explore priorities",
      icon = "\U00002699"))
    q <- length(s$a$Yv)
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "What ω does"),
        tags$p("ω, pronounced omega, is a priority weight. Increasing a variable's weight tells the optimiser to protect precision for that variable, often at some cost to the others.")),
      div(class = "mini-note",
        if (isTRUE(s$density))
          "Density mode re-optimises the analytic design for each weight setting. "
        else "",
        "Re-optimises the full joint design over a grid of priority weights,",
        " trading off two study variables (others keep their relative",
        " weights). Uses a reduced number of multi-starts for speed, so results are",
        " indicative, not publication-grade."),
      fluidRow(
        column(3, selectInput("ox1", "Variable A", s$a$Yv, selected = s$a$Yv[1])),
        column(3, selectInput("ox2", "Variable B", s$a$Yv,
                              selected = s$a$Yv[min(2, q)])),
        column(2, numericInput("ogrid", "Grid points", 7, min = 3, max = 15)),
        column(2, numericInput("ostarts", "Starts", 5, min = 1, max = 25)),
        column(2, actionButton("run_omega", "Sweep",
                               class = "btn-run", style = "margin-top:1.75rem;"))),
      if (has_plotly) plotly_output("pl_omega", height = "430px"),
      DTOutput("dt_omega"),
      div(class = "dl-row",
          downloadButton("dl_omega", "\u03c9 sweep (CSV)", class = "btn-sm"))
    )
  })

  omega_res <- reactiveVal(NULL)

  observeEvent(input$run_omega, {
    s <- strata_res(); req(s)
    if (!require_current_result()) return()
    k1 <- match(input$ox1, s$a$Yv); k2 <- match(input$ox2, s$a$Yv)
    if (is.na(k1) || is.na(k2) || k1 == k2) {
      showNotification("Pick two different study variables.", type = "warning")
      return()
    }
    grid <- seq(0.05, 0.95, length.out = as.integer(input$ogrid))
    total_pair <- s$a$om[k1] + s$a$om[k2]
    out <- vector("list", length(grid))
    withProgress(message = "Sweeping ω…", value = 0, {
      for (i in seq_along(grid)) {
        om <- s$a$om
        om[k1] <- grid[i] * total_pair
        om[k2] <- (1 - grid[i]) * total_pair
        om <- om / sum(om)
        incProgress(1 / length(grid),
                    detail = sprintf("ω(%s) = %.2f", input$ox1, om[k1]))
        if (isTRUE(s$density)) {
          d <- tryCatch(
            omega_sweep_density(
              distributions = s$a$dists, H_vec = s$res$.meta$H_vec,
              n_total = s$settings$n_total, N = s$settings$N,
              ranges = s$a$ranges, omega = s$a$om, pair = c(k1, k2),
              grid = grid[i], n_starts = as.integer(input$ostarts),
              max_iter = 60L,
              constraint = s$settings$constraint,
              cost_total = if (s$settings$constraint == "cost") s$settings$cost_total else NULL,
              cost_fixed = s$settings$cost_fixed,
              cost_per_unit = s$settings$cost_unit,
              cv_targets = s$settings$cv_targets,
              n_min = as.integer(s$settings$n_min)),
            error = function(e) NULL)
          if (!is.null(d)) out[[i]] <- d
        } else {
          d <- tryCatch(
            omega_sweep(
              data = s$a$dat, study_vars = s$a$sv,
              H_vec = s$res$.meta$H_vec, n_total = s$settings$n_total,
              omega = s$a$om, pair = c(k1, k2),
              grid = grid[i], n_starts = as.integer(input$ostarts),
              max_iter = 60L,
              method = s$res$.meta$method,
              constraint = s$settings$constraint,
              cost_total = if (s$settings$constraint == "cost") s$settings$cost_total else NULL,
              cost_fixed = s$settings$cost_fixed,
              cost_per_unit = s$settings$cost_unit,
              cv_targets = s$settings$cv_targets,
              n_min = as.integer(s$settings$n_min),
              distributions = s$res$.meta$distributions),
            error = function(e) NULL)
          if (!is.null(d)) out[[i]] <- d
        }
      }
    })
    omega_res(do.call(rbind, out))
  })

  output$pl_omega <- render_plotly({
    d <- omega_res(); req(d, has_plotly)
    s <- strata_res()
    plot_ly(d, x = ~w) |>
      add_trace(y = ~phi, name = "composite Φ", type = "scatter",
                mode = "lines+markers",
                line = list(color = "#0D1B2A", width = 3),
                marker = list(color = "#0D1B2A")) |>
      add_trace(y = ~cv1, name = paste0("CV ", input$ox1, " (%)"),
                yaxis = "y2", type = "scatter", mode = "lines+markers",
                line = list(color = "#24527A", dash = "dot")) |>
      add_trace(y = ~cv2, name = paste0("CV ", input$ox2, " (%)"),
                yaxis = "y2", type = "scatter", mode = "lines+markers",
                line = list(color = "#2A9D8F", dash = "dot")) |>
      layout(
        xaxis = list(title = paste0("weight on ", input$ox1,
                                    "  (share of the pair's total)")),
        yaxis = list(title = "composite Φ"),
        yaxis2 = list(overlaying = "y", side = "right",
                      title = "CV (%)", showgrid = FALSE),
        legend = list(orientation = "h", y = -0.22),
        margin = list(t = 20, r = 70))
  })

  output$dt_omega <- renderDT({
    d <- omega_res(); req(d)
    names(d) <- c("ω share (A)", "Φ", paste0("CV ", input$ox1, " %"),
                  paste0("CV ", input$ox2, " %"), "gap %")
    datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 15)) |>
      formatRound(1:5, 4)
  })

  output$dl_omega <- downloadHandler(
    filename = function() "multistrata_omega_sweep.csv",
    content = function(f) {
      d <- omega_res(); req(d)
      utils::write.csv(d, f, row.names = FALSE)
    })

  # ============================================================================
  # Deploy tab: take the design into production
  #
  # The optimisation tabs stop at a set of boundaries and an allocation. This
  # tab covers the step after that, using the package's own operational
  # helpers: assign_strata() scores an arbitrary frame into the fitted strata,
  # and as_svydesign() turns a drawn sample into a design object ready for
  # design-based estimation with survey or srvyr.
  # ============================================================================
  output$ui_deploy <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      paste("Run a design first, then use this tab to score a frame into the",
            "optimal strata and export a ready-to-analyse survey design."),
      title = "Run first to deploy",
      icon = "\U0001F69B"))
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "What this does"),
        tags$p(paste(
          "A design is only useful once it is applied to units. This tab",
          "labels every unit of a frame with the stratum it falls in, draws a",
          "stratified sample using the optimal allocation, and exports the",
          "result with the design weights already attached, so it can go",
          "straight into estimation."))),
      div(class = "mini-note",
        tags$strong("assign_strata()"),
        " scores the frame; ", tags$strong("as_svydesign()"),
        " attaches .stratum, .N_h, .n_h and .weight (= Nₕ/nₕ) to the drawn",
        " sample. Both are package functions, so the downloads below are",
        " reproducible from the R Code tab."),
      fluidRow(
        column(5,
          radioButtons("dep_src", "Frame to score",
            choices = c("Fitting population" = "self", "Upload a frame" = "up"),
            selected = if (isTRUE(s$density)) "up" else "self")),
        column(7,
          conditionalPanel("input.dep_src == 'up'",
            fileInput("dep_file",
                      if (has_readxl) "CSV or Excel frame" else "CSV frame",
                      accept = if (has_readxl)
                        c(".csv", ".xlsx", ".xls") else ".csv"),
            div(class = "mini-note",
                "Must contain the stratification variables ",
                tags$strong(paste(s$a$Xv, collapse = ", ")),
                " (matched by name), and the study variables too if you want",
                " a design object rather than just stratum labels.")))),
      if (isTRUE(s$density))
        div(class = "mini-note",
            "This is a density-mode design, so there is no fitting",
            " population: upload the frame you intend to survey."),
      fluidRow(
        column(4, numericInput("dep_seed", "Sampling seed", value = 1,
                               min = 1, step = 1)),
        column(4, checkboxInput("dep_planned",
                                "Weight by planned nₕ", value = FALSE)),
        column(4, actionButton("run_deploy", "\u25B6  Score and draw",
                               class = "btn-run",
                               style = "width:auto; margin-top:1.7rem;"))),
      uiOutput("ui_deploy_out"))
  })

  deploy_res <- reactiveVal(NULL)
  observeEvent(input$dep_src, deploy_res(NULL))
  observeEvent(input$dep_file, deploy_res(NULL))

  observeEvent(input$run_deploy, {
    s <- strata_res(); req(s)
    if (!require_current_result()) return()

    frame <- tryCatch({
      if (identical(input$dep_src, "up")) {
        if (is.null(input$dep_file))
          stop("Choose a frame file first.")
        ext <- tolower(tools::file_ext(input$dep_file$name))
        if (ext == "csv") utils::read.csv(input$dep_file$datapath)
        else if (ext %in% c("xlsx", "xls")) {
          if (!has_readxl)
            stop("Install the 'readxl' package to read Excel files.")
          as.data.frame(readxl::read_excel(input$dep_file$datapath))
        } else stop("Unsupported file type: .", ext)
      } else {
        if (isTRUE(s$density))
          stop("A density-mode design has no fitting population: upload a frame.")
        as.data.frame(dataset())
      }
    }, error = function(e) e)

    if (inherits(frame, "error")) {
      deploy_res(list(err = conditionMessage(frame)))
      return()
    }

    out <- tryCatch({
      withProgress(message = "Scoring the frame\u2026", value = 0.3, {
        # assign_strata() matches the stratification columns by name and
        # errors informatively when they are absent, so let it validate.
        sid <- assign_strata(s$res, frame, labels = FALSE)
        lab <- assign_strata(s$res, frame, labels = TRUE)
        incProgress(0.3, detail = "drawing the stratified sample")

        alloc <- s$res$allocation
        set.seed(as.integer(input$dep_seed %||% 1))
        take <- unlist(lapply(seq_along(alloc), function(h) {
          idx <- which(sid == h)
          if (!length(idx)) return(integer(0))
          sample(idx, min(length(idx), as.integer(round(alloc[h]))))
        }))

        scored <- frame
        scored$.stratum <- lab
        incProgress(0.3, detail = "attaching design weights")

        # frame_only avoids a hard dependency on survey for the download;
        # the svydesign object itself is built only for the on-screen summary.
        wframe <- as_svydesign(s$res, frame[take, , drop = FALSE],
                               use_planned = isTRUE(input$dep_planned),
                               frame_only = TRUE)
        des_txt <- if (requireNamespace("survey", quietly = TRUE))
          tryCatch(paste(utils::capture.output(print(
            as_svydesign(s$res, frame[take, , drop = FALSE],
                         use_planned = isTRUE(input$dep_planned)))),
            collapse = "\n"), error = function(e) NULL) else NULL

        ss <- s$res$stratum_stats
        cover <- data.frame(
          stratum = ss$stratum,
          N_h_design = ss$N_h,
          n_h_planned = as.integer(round(alloc)),
          N_h_frame = as.integer(table(factor(lab, levels = ss$stratum))),
          n_h_drawn = as.integer(table(factor(lab[take], levels = ss$stratum))),
          stringsAsFactors = FALSE)
        cover$shortfall <- cover$n_h_planned - cover$n_h_drawn

        list(scored = scored, wframe = wframe, cover = cover,
             des_txt = des_txt, n_frame = nrow(frame), n_drawn = length(take))
      })
    }, error = function(e) list(err = conditionMessage(e)))

    deploy_res(out)
  })

  output$ui_deploy_out <- renderUI({
    r <- deploy_res(); if (is.null(r)) return(NULL)
    if (!is.null(r$err))
      return(div(class = "err-box", "\u2716 ", r$err))

    short <- sum(pmax(0L, r$cover$shortfall))
    tagList(
      div(class = if (short > 0) "warn-box" else "ok-box",
        tags$strong(sprintf("Scored %s units; drew %s.",
                            format(r$n_frame, big.mark = ","),
                            format(r$n_drawn, big.mark = ","))),
        if (short > 0) sprintf(
          paste("  %d sample units could not be drawn because those strata",
                "hold fewer frame units than the allocation asks for; see",
                "the shortfall column."), short) else
          "  Every stratum supplied its full allocation."),
      tags$h6(style = "margin-top:1rem;", "Stratum coverage"),
      div(class = "mini-note",
          "Design Nₕ comes from the fitting population; frame Nₕ from the",
          " frame you just scored. Large differences mean the frame's",
          " composition has drifted from the population the design was",
          " optimised for."),
      renderDT(datatable(r$cover, rownames = FALSE,
                         options = list(dom = "t", pageLength = 50,
                                        scrollX = TRUE))),
      if (!is.null(r$des_txt)) tagList(
        tags$h6(style = "margin-top:1rem;", "survey design object"),
        div(class = "code-box", tags$pre(r$des_txt))),
      tags$h6(style = "margin-top:1rem;", "Weighted sample (first rows)"),
      renderDT(datatable(utils::head(r$wframe, 20), rownames = FALSE,
                         options = list(dom = "t", scrollX = TRUE))),
      div(class = "dl-row",
        downloadButton("dl_scored", "Scored frame (CSV)", class = "btn-sm"),
        downloadButton("dl_wframe", "Weighted sample (CSV)", class = "btn-sm"),
        downloadButton("dl_cover", "Stratum coverage (CSV)", class = "btn-sm")))
  })

  output$dl_scored <- downloadHandler(
    filename = function() "multistrata_scored_frame.csv",
    content = function(f) {
      r <- deploy_res(); req(r, is.null(r$err))
      utils::write.csv(r$scored, f, row.names = FALSE)
    })
  output$dl_wframe <- downloadHandler(
    filename = function() "multistrata_weighted_sample.csv",
    content = function(f) {
      r <- deploy_res(); req(r, is.null(r$err))
      utils::write.csv(r$wframe, f, row.names = FALSE)
    })
  output$dl_cover <- downloadHandler(
    filename = function() "multistrata_stratum_coverage.csv",
    content = function(f) {
      r <- deploy_res(); req(r, is.null(r$err))
      utils::write.csv(r$cover, f, row.names = FALSE)
    })

  # ============================================================================
  # Sample size tab (find_n): smallest n meeting CV targets
  # ============================================================================
  output$ui_findn <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run a design first, then find the smallest sample size that meets your CV precision targets.",
      title = "Run first to search sample size",
      icon = "\U0001F50D"))
    q <- length(s$a$Yv)
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "What this answers"),
        tags$p("This tool asks: what is the smallest total sample size that still reaches my target precision? It reuses the current strata and searches for the smallest n that satisfies your CV targets.")),
      div(class = "mini-note",
        "The smallest total sample size n that keeps every study variable's ",
        "coefficient of variation (CV) at or below your target, reusing the ",
        "current strata. Enter a target percentage per variable and press Find."),
      fluidRow(lapply(seq_len(q), function(k)
        column(max(3, floor(12 / q)),
          numericInput(paste0("fn_cv_", k),
                       paste0("CV target · ", s$a$Yv[k], " (%)"),
                       value = 5, min = 0.1, max = 100, step = 0.5)))),
      actionButton("run_findn", "▶  Find smallest n", class = "btn-run",
                   style = "width:auto; margin-top:0.3rem;"),
      uiOutput("ui_findn_out"))
  })

  findn_res <- reactiveVal(NULL)
  observeEvent(input$run_findn, {
    s <- strata_res(); req(s)
    if (!require_current_result()) return()
    q   <- length(s$a$Yv)
    cvt <- vapply(seq_len(q),
                  function(k) as.numeric(input[[paste0("fn_cv_", k)]]) / 100, 1)
    withProgress(message = "Searching for the smallest n…", value = 0.4, {
      fn <- tryCatch(
        if (isTRUE(s$density))
          find_n_density(result = s$res, cv_targets = cvt)
        else
          find_n(s$a$dat, s$a$sv, cv_targets = cvt, result = s$res),
        error = function(e) e)
      findn_res(list(fn = fn, cvt = cvt))
    })
  })

  output$ui_findn_out <- renderUI({
    r <- findn_res(); if (is.null(r)) return(NULL)
    if (inherits(r$fn, "error"))
      return(div(class = "err-box", "✖ ", conditionMessage(r$fn)))
    fn <- r$fn; s <- strata_res()
    if (!isTRUE(fn$feasible))
      return(div(class = "warn-box", tags$strong("Targets not attainable "),
        "even at a full census with these strata. Try looser CV targets, or ",
        "add strata on the 'How many strata' tab."))
    cv_txt <- paste(sprintf("%s: %.2f%% (target %.2f%%)", s$a$Yv,
                            100 * fn$cv_achieved, 100 * r$cvt), collapse = "  ·  ")
    tagList(
      div(class = "ok-box",
        tags$strong(sprintf("Smallest feasible n = %d", fn$n)),
        sprintf("   (searched %d to %d; %d designs evaluated)",
                fn$n_range[1], fn$n_range[2], nrow(fn$trace))),
      div(class = "bounds-card", tags$strong("Achieved CV:  "), cv_txt),
      div(class = "mini-note",
        sprintf("Composite variance at this n: %s.  The allocation is shown on ",
                fmt_num(fn$composite_variance, 8)),
        "the Optimal Design tab if you re-run with n = ", tags$strong(fn$n), "."),
      div(class = "dl-row",
          downloadButton("dl_findn_trace", "Sample-size trace (CSV)", class = "btn-sm"),
          downloadButton("dl_findn_alloc", "Feasible allocation (CSV)", class = "btn-sm")))
  })

  output$dl_findn_trace <- downloadHandler(
    filename = function() "multistrata_sample_size_trace.csv",
    content = function(f) {
      r <- findn_res(); req(r, inherits(r$fn, "findN"))
      utils::write.csv(r$fn$trace, f, row.names = FALSE)
    })

  output$dl_findn_alloc <- downloadHandler(
    filename = function() "multistrata_sample_size_allocation.csv",
    content = function(f) {
      r <- findn_res(); req(r, inherits(r$fn, "findN"))
      fn <- r$fn
      out <- data.frame(stratum = seq_along(fn$allocation),
                        n_h = as.integer(fn$allocation),
                        stringsAsFactors = FALSE)
      utils::write.csv(out, f, row.names = FALSE)
    })

  # ============================================================================
  # How many strata tab (choose_H): trace Phi* across candidate grids
  # ============================================================================
  output$ui_chooseH <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run a design first, then explore how precision improves as you add more strata.",
      title = "Run first to explore strata counts",
      icon = "\U0001F5C2"))
    tagList(
      div(class = "explain-card",
        div(class = "explain-title", "What this answers"),
        tags$p("This tool asks: would more strata noticeably improve precision? The recommended point is an elbow: after that, extra strata give smaller gains and may be harder to manage in the field.")),
      div(class = "mini-note",
        if (isTRUE(s$density))
          "Density mode traces the analytic optimum for each candidate grid. "
        else "",
        "How the optimal composite variance Φ* falls as you use more strata per ",
        "variable. The marked point is the elbow. Beyond it, extra strata buy ",
        "little precision for the added field cost."),
      fluidRow(
        column(5, sliderInput("ch_range", "Strata per variable to try",
                              min = 2, max = 8, value = c(2, 5), step = 1)),
        column(3, numericInput("ch_starts", "Multi-starts", 6, min = 1, max = 25)),
        column(4, actionButton("run_chooseH", "▶  Trace Φ*", class = "btn-run",
                               style = "width:auto; margin-top:1.7rem;"))),
      plotOutput("pl_chooseH", height = "380px"),
      DTOutput("dt_chooseH"),
      uiOutput("ui_chooseH_rec"),
      div(class = "dl-row",
          downloadButton("dl_chooseH", "Strata trace (CSV)", class = "btn-sm")))
  })

  chooseH_res <- reactiveVal(NULL)
  observeEvent(input$run_chooseH, {
    s <- strata_res(); req(s)
    if (!require_current_result()) return()
    rng <- input$ch_range
    withProgress(message = "Tracing Φ* across grids…", value = 0.3, {
      ch <- tryCatch(
        if (isTRUE(s$density))
          choose_H_density(
            distributions = s$a$dists, n_total = as.integer(s$settings$n_total),
            N = as.integer(s$settings$N), H_each = seq.int(rng[1], rng[2]),
            ranges = s$a$ranges, omega = s$a$om,
            constraint = s$settings$constraint,
            cost_total = if (s$settings$constraint == "cost") s$settings$cost_total else NULL,
            cost_fixed = s$settings$cost_fixed,
            cost_per_unit = s$settings$cost_unit,
            cv_targets = s$settings$cv_targets,
            n_min = s$settings$n_min,
            n_starts = as.integer(input$ch_starts), plot = FALSE)
        else
          choose_H(s$a$dat, s$a$sv, n_total = as.integer(s$settings$n_total),
                   H_each = seq.int(rng[1], rng[2]), omega = s$a$om,
                   distributions = s$a$dists,
                   constraint = s$settings$constraint,
                   cost_total = if (s$settings$constraint == "cost") s$settings$cost_total else NULL,
                   cost_fixed = s$settings$cost_fixed,
                   cost_per_unit = s$settings$cost_unit,
                   cv_targets = s$settings$cv_targets,
                   n_min = s$settings$n_min,
                   method = s$settings$method,
                   n_starts = as.integer(input$ch_starts), plot = FALSE),
        error = function(e) e)
      chooseH_res(ch)
    })
  })

  output$pl_chooseH <- renderPlot({
    ch <- chooseH_res(); req(ch); if (inherits(ch, "error")) return(NULL)
    tryCatch(
      plot(ch),
      error = function(e) {
        op <- graphics::par(mar = c(1, 1, 1, 1))
        on.exit(graphics::par(op), add = TRUE)
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Strata trace plot unavailable in this view.")
      })
  })
  output$dt_chooseH <- renderDT({
    ch <- chooseH_res(); req(ch); if (inherits(ch, "error")) return(NULL)
    tb <- ch$table
    names(tb) <- c("Grid", "Total strata", "Φ*", "Optimality gap",
                   "RE vs coarsest", "Marginal gain")
    datatable(tb, rownames = FALSE, options = list(dom = "t")) |>
      formatRound(c("Φ*"), 6) |>
      formatPercentage(c("Optimality gap", "Marginal gain"), 1) |>
      formatRound("RE vs coarsest", 3)
  })
  output$ui_chooseH_rec <- renderUI({
    ch <- chooseH_res()
    if (is.null(ch) || inherits(ch, "error")) {
      if (inherits(ch, "error")) return(div(class = "err-box", conditionMessage(ch)))
      return(NULL)
    }
    div(class = "ok-box", tags$strong("Recommended: "),
        sprintf("%s strata (%d in total), the elbow of the curve.",
                paste(ch$recommended, collapse = " × "), prod(ch$recommended)))
  })

  output$dl_chooseH <- downloadHandler(
    filename = function() "multistrata_strata_trace.csv",
    content = function(f) {
      ch <- chooseH_res(); req(ch, inherits(ch, "chooseH"))
      utils::write.csv(ch$table, f, row.names = FALSE)
    })

  # ============================================================================
  # Report tab: shareable self-contained HTML summary of the design
  # ============================================================================
  output$ui_report <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run a design first, then download a shareable report for colleagues or supervisors.",
      title = "Run first to create a report",
      icon = "\U0001F4C4"))
    tagList(
      div(class = "action-card",
        div(class = "action-card-title", "Shareable design report"),
        div(class = "mini-note",
          "Download a self-contained HTML report with the design settings,",
          " boundaries, allocation, precision results, figures and any",
          " follow-on analyses you have run."),
        div(class = "dl-row",
          downloadButton("dl_report", "⬇  Download HTML report",
                         class = "btn-run", style = "width:auto;"),
          if (has_pagedown)
            downloadButton("dl_report_pdf", "⬇  Download PDF report",
                           class = "btn-run", style = "width:auto;"))),
      div(class = "mini-note",
        if (has_pagedown)
          "Tip: the HTML report also includes a Save as PDF button."
        else paste("Tip: open the HTML report and press its ‘Save as PDF’ button",
                   "(or your browser’s Print → Save as PDF) to get a PDF. For a",
                   "direct PDF button here, install the ‘pagedown’ package.")))
  })

  report_state <- reactive({
    s <- strata_res(); req(s)
    out <- s
    out$analysis <- list(
      benchmark = tryCatch(bench_res(), error = function(e) NULL),
      findn     = findn_res(),
      chooseH   = {
        ch <- chooseH_res()
        if (inherits(ch, "error")) NULL else ch
      },
      omega     = omega_res()
    )
    out
  })

  output$dl_report <- downloadHandler(
    filename = function() "multistrata_design_report.html",
    content  = function(f) writeLines(.build_report_html(report_state()), f))

  output$dl_report_pdf <- downloadHandler(
    filename = function() "multistrata_design_report.pdf",
    content  = function(f) {
      tmp <- tempfile(fileext = ".html")
      writeLines(.build_report_html(report_state()), tmp)
      tryCatch(
        pagedown::chrome_print(tmp, output = f),
        error = function(e)
          stop("PDF export needs Chrome/Chromium (via the 'pagedown' package). ",
               "Use the HTML report's 'Save as PDF' button instead. Details: ",
               conditionMessage(e)))
    })

  # ============================================================================
  # TAB 7: R Code
  # ============================================================================
  gen_code <- reactive({
    s <- strata_res(); req(s)
    a <- s$a; st <- s$settings

    # Density (data-free) mode emits a stratify_multi_density() script.
    if (isTRUE(s$density)) {
      dist_lines <- vapply(seq_along(a$dists), function(j) {
        d <- a$dists[[j]]; pars <- d[setdiff(names(d), "family")]
        sprintf('  %s = list(family = "%s", %s)', a$Xv[j], d$family,
                paste(sprintf("%s = %s", names(pars),
                              vapply(pars, format, "")), collapse = ", "))
      }, "")
      range_lines <- vapply(seq_along(a$ranges), function(j)
        sprintf("c(%s, %s)", format(a$ranges[[j]][1]),
                format(a$ranges[[j]][2])), "")
      cost_block <- if (st$constraint == "cost")
        sprintf(paste0(',\n  constraint    = "cost",',
                       '\n  cost_total    = %s, cost_fixed = %s,',
                       '\n  cost_per_unit = %s'),
                format(st$cost_total), format(st$cost_fixed),
                format(st$cost_unit)) else ""
      return(paste0(
        "# Reproducible script (density mode), multistrata app\n",
        "# ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n",
        "library(multistrata)\n\n",
        "# The multi-start solver draws random initial boundaries, so this\n",
        "# seed is what makes the design below identical to the one the app\n",
        "# displayed. Change it to explore restart sensitivity.\n",
        sprintf("set.seed(%d)\n\n", as.integer(st$seed %||% 1L)),
        "res <- stratify_multi_density(\n",
        "  distributions = list(\n", paste(dist_lines, collapse = ",\n"), "\n  ),\n",
        sprintf("  H_vec         = c(%s),\n", paste(a$Hs, collapse = ", ")),
        sprintf("  n_total       = %d,\n", st$n_total),
        sprintf("  N             = %d,\n", st$N),
        sprintf("  ranges        = list(%s),\n", paste(range_lines, collapse = ", ")),
        sprintf("  omega         = c(%s),\n", paste(round(a$om, 4), collapse = ", ")),
        sprintf("  n_min         = %d,\n", st$n_min),
        sprintf("  n_starts      = %d,\n", st$n_starts),
        sprintf("  max_iter      = %d", st$max_iter),
        cost_block,
        "\n)\n\nprint(res)\nsummary(res)\n\n",
        "# Density benchmark: joint optimum vs simpler analytic alternatives\n",
        "bench <- benchmark_density_designs(result_joint = res)\n",
        "print(bench)\n\n",
        "# Precision-driven sample size example: edit targets as needed\n",
        sprintf("fn <- find_n_density(result = res, cv_targets = c(%s))\n",
                paste(rep("0.05", length(a$Yv)), collapse = ", ")),
        "print(fn)\n\n",
        "# Number-of-strata trace example\n",
        "ch <- choose_H_density(\n",
        "  distributions = res$.meta$distributions,\n",
        "  n_total       = res$.meta$n_total,\n",
        "  N             = res$.meta$N,\n",
        "  H_each        = 2:5,\n",
        "  ranges        = res$ranges,\n",
        "  omega         = res$.meta$omega,\n",
        "  n_starts      = 5,\n",
        "  plot          = FALSE\n",
        ")\n",
        "print(ch)\n\n",
        "# Omega explorer example for the first two variables\n",
        "om <- omega_sweep_density(\n",
        "  distributions = res$.meta$distributions,\n",
        "  H_vec         = res$.meta$H_vec,\n",
        "  n_total       = res$.meta$n_total,\n",
        "  N             = res$.meta$N,\n",
        "  ranges        = res$ranges,\n",
        "  omega         = res$.meta$omega,\n",
        "  pair          = c(1, 2),\n",
        "  n_starts      = 5\n",
        ")\n",
        "print(om)\n"))
    }

    data_block <- switch(st$data_source,
      builtin = sprintf(
        'library(multistrata)\n\ndata(%s)\ndf <- %s',
        builtin_specs[[st$builtin_ds]]$ds,
        builtin_specs[[st$builtin_ds]]$ds),
      upload = sprintf(
        'library(multistrata)\n\ndf <- read.csv("%s")   # adjust path',
        st$file_name %||% "your_data.csv"),
      simulate = sprintf(
        'library(multistrata)\n\ndf <- sim_population(N = %d, p = %d, q = %d,\n                     scenario = "%s", seed = %d)',
        st$sim_N, st$sim_p, st$sim_q, st$sim_scenario, st$sim_seed))

    dist_block <- ""
    if (isTRUE(a$auto)) {
      # emitted inline as distributions = "auto"
    } else if (!is.null(a$dists)) {
      lines <- vapply(seq_along(a$dists), function(j) {
        d <- a$dists[[j]]
        pars <- d[setdiff(names(d), "family")]
        sprintf('  %s = list(family = "%s", %s)', a$Xv[j], d$family,
                paste(sprintf("%s = %s", names(pars),
                              vapply(pars, format, "")), collapse = ", "))
      }, "")
      dist_block <- sprintf("\ndists <- list(\n%s\n)\n",
                            paste(lines, collapse = ",\n"))
    }

    cv_block <- if (!is.null(a$cvt))
      sprintf(",\n  cv_targets    = c(%s)",
              paste(format(a$cvt), collapse = ", ")) else ""
    cost_block <- if (st$constraint == "cost")
      sprintf(paste0(',\n  constraint    = "cost",',
                     '\n  cost_total    = %s, cost_fixed = %s,',
                     '\n  cost_per_unit = %s'),
              format(st$cost_total), format(st$cost_fixed),
              format(st$cost_unit)) else ""

    # Deployment step: score a frame into the fitted strata and export a
    # design object. plot_strata_rgl() is emitted rather than rendered,
    # because it needs an interactive OpenGL device that Shiny cannot host.
    deploy_block <- paste0(
      "# ---------------------------------------------------------------\n",
      "# Deploy the design\n",
      "# ---------------------------------------------------------------\n",
      "# 1. Score a frame into the optimal strata. The fitting population is\n",
      "#    reused here; swap in a fresh census extract or hold-out frame.\n",
      "scored <- assign_strata(res, X)\n",
      "table(scored)\n\n",
      "# 2. Draw a stratified sample using the optimal allocation.\n",
      "set.seed(1)\n",
      "sid  <- assign_strata(res, X, labels = FALSE)\n",
      "take <- unlist(lapply(seq_along(res$allocation), function(h) {\n",
      "  idx <- which(sid == h)\n",
      "  sample(idx, min(length(idx), res$allocation[h]))\n",
      "}))\n\n",
      "# 3. Export the drawn sample as a ready-to-analyse design object.\n",
      "#    as_svydesign() expects the SAMPLED units, not the population.\n",
      "if (requireNamespace(\"survey\", quietly = TRUE)) {\n",
      "  des <- as_svydesign(res, df[take, ])\n",
      "  print(des)\n",
      "}\n\n",
      "# ... or the augmented frame alone, carrying .stratum, .N_h, .n_h and\n",
      "# .weight, for srvyr and other tools.\n",
      "frame <- as_svydesign(res, df[take, ], frame_only = TRUE)\n",
      "head(frame)\n",
      if (length(a$Xv) == 2L) paste0(
        "\n# Interactive 3-D boundary view: needs the 'rgl' package and an OpenGL\n",
        "# device, so it runs in R rather than inside the app.\n",
        "if (requireNamespace(\"rgl\", quietly = TRUE))\n",
        "  plot_strata_rgl(res, data = X, study_vars = Y)\n") else "")

    paste0(
      "# Reproducible script generated by the multistrata app\n",
      "# ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n",
      data_block, "\n\n",
      "X <- df[, c(", paste(sprintf('"%s"', a$Xv), collapse = ", "), ")]\n",
      "Y <- df[, c(", paste(sprintf('"%s"', a$Yv), collapse = ", "), ")]\n",
      dist_block,
      "\n# The multi-start solver draws random initial boundaries, so this seed\n",
      "# is what makes the design below identical to the one the app displayed.\n",
      "# Change it to explore restart sensitivity.\n",
      sprintf("set.seed(%d)\n", as.integer(st$seed %||% 1L)),
      "\nres <- stratify_multi(\n",
      "  data          = X,\n",
      "  study_vars    = Y,\n",
      sprintf("  H_vec         = c(%s),\n",
              paste(s$res$.meta$H_vec, collapse = ", ")),
      sprintf("  n_total       = %d,\n", st$n_total),
      sprintf("  omega         = c(%s),\n",
              paste(round(a$om, 4), collapse = ", ")),
      sprintf("  n_min         = %d,\n", st$n_min),
      sprintf("  n_starts      = %d,\n", st$n_starts),
      sprintf("  max_iter      = %d,\n", st$max_iter),
      sprintf('  method        = "%s"', st$method),
      if (isTRUE(a$auto)) ',\n  distributions = "auto"'
      else if (!is.null(a$dists)) ",\n  distributions = dists" else "",
      cost_block, cv_block,
      "\n)\n\nprint(res)\nsummary(res)\n\n",
      "# Design figures\n",
      "plot_strata(res, data = X, study_vars = Y)\n",
      "plot(res, which = \"allocation\")\n\n",
      "# Benchmark against sequential / proportional / equal designs\n",
      sprintf(paste0("bench <- benchmark_designs(X, Y, H_vec = c(%s),\n",
                     "                           n_total = %d, omega = c(%s),\n",
                     "                           result_joint = res)\n"),
              paste(s$res$.meta$H_vec, collapse = ", "), sum(s$res$allocation),
              paste(round(a$om, 4), collapse = ", ")),
      "print(bench)\n\n",
      "# Precision-driven sample size example: edit targets as needed\n",
      sprintf("fn <- find_n(X, Y, result = res, cv_targets = c(%s))\n",
              paste(rep("0.05", length(a$Yv)), collapse = ", ")),
      "print(fn)\n\n",
      "# Number-of-strata trace example\n",
      "ch <- choose_H(X, Y,\n",
      "  n_total  = res$.meta$n_total,\n",
      "  H_each   = 2:5,\n",
      "  omega    = res$.meta$omega,\n",
      "  n_starts = 5,\n",
      "  plot     = FALSE\n",
      ")\n",
      "print(ch)\n\n",
      sprintf(paste0("# Omega explorer: precision trade-off between %s and %s\n"),
              a$Yv[1], a$Yv[min(2L, length(a$Yv))]),
      "om <- omega_sweep(X, Y,\n",
      "  n_total  = res$.meta$n_total,\n",
      "  result   = res,\n",
      sprintf("  pair     = c(%s),\n",
              paste(sprintf('"%s"', a$Yv[seq_len(min(2L, length(a$Yv)))]),
                    collapse = ", ")),
      "  n_starts = 5\n",
      ")\n",
      "print(om)\n\n",
      deploy_block)
  })

  output$ui_code <- renderUI({
    s <- strata_res(); if (is.null(s)) return(empty_state(
      "Run the optimisation first, then this tab shows a reproducible R script for advanced users.",
      title = "Run first to create code",
      icon = "\U0001F4BB"))
    tagList(
      div(class = "action-card",
        div(class = "action-card-title", "Reproducible R script"),
        div(class = "mini-note",
            "Copy or download this script to reproduce the current design",
            " with the multistrata package alone."),
        downloadButton("dl_code", "⬇  Download R script",
                       class = "btn-run", style = "width:auto;")),
      div(class = "code-box", verbatimTextOutput("code_txt")))
  })
  output$code_txt <- renderText(gen_code())
  output$dl_code <- downloadHandler(
    filename = function() "multistrata_design.R",
    content = function(f) writeLines(gen_code(), f))
}

shinyApp(ui, server)
