# multistrata Shiny application

This interactive application ships **inside** the `multistrata` package
(installed under `inst/shiny/`). It is the companion to the methodology and to
the package's exported functions, everything it computes is done by the
package itself.

## Launching

From any R session with the package installed:

```r
library(multistrata)
run_app()                       # standalone window, or the default browser
run_app(launch.browser = TRUE)  # force the default web browser
run_app(port = 8123)            # fixed port
```

`run_app()` requires the interface packages `shiny`, `bslib`, `DT` and
`jsonlite` (declared in the package's `Suggests`). Optional packages add
features: `plotly` (interactive plots), `readxl` (Excel upload), `survey`
(design objects on the Deploy tab), `pagedown` (direct PDF report export) and
`FitVerse` (an extra distribution-fitting engine). Each is used conditionally:
without `plotly`, for example, the Visualise tab falls back to the static
package plots rather than failing. Install the required set once with:

```r
install.packages(c("shiny", "bslib", "DT", "jsonlite"))
```

## Tabs

| Tab | Package functions behind it |
|---|---|
| Introduction | — |
| Optimal Design | `stratify_multi()`, `stratify_multi_density()`, `fit_distributions()` |
| Visualise | `plot_strata()`, `plot()`, `summary()` |
| Benchmark | `benchmark_designs()`, `benchmark_density_designs()`, `relative_efficiency()` |
| Precision (CV) | `design_cv()`, `composite_variance()`, `variance_lower_bound()` |
| Sample size | `find_n()`, `find_n_density()` |
| No. of Strata | `choose_H()`, `choose_H_density()` |
| ω Explorer | `omega_sweep()`, `omega_sweep_density()` |
| Deploy | `assign_strata()`, `as_svydesign()` |
| Report | — (renders the above) |
| R Code | emits a runnable script using the functions above |

The Deploy tab covers the step after optimisation: it scores a frame (the
fitting population, or an uploaded CSV/Excel frame) into the fitted strata
with `assign_strata()`, draws a stratified sample using the optimal
allocation, and exports it through `as_svydesign()` with `.stratum`, `.N_h`,
`.n_h` and `.weight` attached, plus a stratum-coverage table flagging any
stratum whose frame holds fewer units than the allocation requires.
`plot_strata_rgl()` needs an interactive OpenGL device, so the R Code tab
emits the call rather than the app rendering it.

Every exported function is reachable from the app except
`compute_stratum_stats_parametric()`, which computes analytic stratum moments
at a given set of boundaries. The parametric mode calls it internally on every
run, but it is a building block for users assembling their own objective
rather than a workflow step, so it has no control of its own.

## Files

- `app.R`, the single-file Shiny app (UI + server).
- `app_helpers.R`, input-validation and snapshot helpers, sourced by `app.R`.
- `www/logo.svg`, the app logo.
- `khan_reddy_case1.csv`, `khan_reddy_case2.csv`, example upload datasets
  reproducing the bivariate Khan & Reddy population.
- `reproduce.R`, optional console utilities (regenerate the example CSVs and
  reproduce the paper's Case 1 / Case 2 designs). Not loaded on launch:
  `source(system.file("shiny", "reproduce.R", package = "multistrata"))`.

## Development

When working from the package source (not an installed copy), load it with
`pkgload::load_all()` and then call `run_app()`, `system.file()` resolves to
`inst/shiny/` under `load_all()`, so the same launcher works in development.
