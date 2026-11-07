# multistrata 1.0-1

First release on CRAN, and the first under this name: the package was
developed as `stratifyMultiR` and renamed to `multistrata` before release
(all lower case, easier to type, and it reads as what it does). Two names
changed with it and nothing else: the package itself, and the class of a
fitted design, which is now `"multistrata"` rather than `"stratifyMulti"`.
Every function and argument name is unchanged, so a script written against
the development version needs only `library(multistrata)`, plus, if it
tested the class explicitly, `inherits(res, "multistrata")`.

The entries below describe the package relative to the pre-release
development versions.

## New features

* `omega_sweep()`, traces the precision trade-off between two study variables
  as the priority weight shifts from one to the other, re-optimising the design
  at each grid point. This is the data-mode counterpart of
  `omega_sweep_density()`, completing the data/density pairing that
  `find_n()`/`find_n_density()`, `choose_H()`/`choose_H_density()` and
  `benchmark_designs()`/`benchmark_density_designs()` already follow. The
  swept pair may be given by name or by index, and `result=` reuses a fitted
  design's `H_vec`, `omega`, `method` and `distributions`. Both sweep
  functions now return the swept indices and names as the `"pair"` and
  `"pair_names"` attributes.

* `stratify_multi_density()`, data-free ("analytic") stratification for
  designs where the stratification variables are themselves the study
  variables. It takes only the marginal
  distributions (no population data frame): stratum weights come from the
  marginal CDFs and stratum means/variances from Gauss–Kronrod quadrature of
  the marginal densities (`stats::integrate`). This is the exact analogue of
  the density-based LINGO formulation and matches the empirical mode of
  `stratify_multi()` in the large-N limit (verified: analytic moments equal
  large-sample empirical moments to 3–4 dp).

* `find_n()`, precision-driven ("inverse") design. Returns the smallest total
  sample size for which an allocation exists that meets a coefficient-of-
  variation target on every study variable, by integer bisection over the
  CV-constrained allocation sub-problem. Boundaries may be reused from a fitted
  object (`result=`) or optimised internally; `reoptimize=TRUE` re-optimises
  boundaries at each trial size.

* `assign_strata()`, scores an arbitrary frame (estimation frame, fresh census
  extract, hold-out data) into the optimal strata of a fitted design. Matches
  stratification columns by name when available.

* `as_svydesign()`, exports a ready-to-analyse `survey::svydesign` object (or,
  with `frame_only=TRUE`, an augmented data frame carrying `.stratum`, `.N_h`,
  `.n_h` and `.weight` for `srvyr` and other tools).

* `choose_H()`, traces the optimal composite variance against a family of
  candidate strata grids and flags the elbow, answering "how many strata?".
  Comes with `print` and `plot` methods.

* `stratify_multi()` gains a `trace` argument; when `TRUE` the returned object
  carries `phi_trace`, the per-iteration objective sequence of the best start
  (used to visualise the monotone descent of Section 4.2).

* `sim_population()` gains two scenarios: `"weibull_pareto"` (Scenario C,
  heavy-tailed establishment-survey archetype) and `"copula"` (Scenario E,
  dependent stratification variables via a Gaussian copula with log-normal
  margins).

## Bug fixes

* `stratify_multi(cv_targets = ...)` now works. The CV-constrained allocation
  wrote its finite-population correction as `1/n - 1/N_h`, which is `Inf` for
  an empty stratum; since such a stratum also has `W_h = 0`, the variance term
  became `0 * Inf = NaN`, and the `NaN` reached the SLSQP solver as its
  warm-start objective and aborted the run with *"objective in x0 returns
  NA"*. Empty strata arise routinely while the boundary search explores trial
  cuts, so the failure became more likely as `n_starts` grew and was
  near-certain at the default of 50, the coefficient-of-variation constraint
  was effectively unusable. The correction is now guarded exactly as
  `composite_variance()`, `design_cv()` and the internal objective already
  guarded it. Reported CVs and allocations are unchanged wherever the function
  previously succeeded.

* `plot_strata_rgl()` now tells the truth when `rgl` is present but will not
  load. Recent macOS releases ship no OpenGL, so `rgl`'s dlopen of `libGLU`
  fails; the old message said *"Package 'rgl' is required ... install it"*,
  which is useless advice to someone who already has it. The two cases are
  now distinguished, and the unloadable case points at
  `plot_strata()`, whose two-variable view draws the same 3-D perspective
  plot with boundary walls and needs no OpenGL. `rgl`'s own dlopen warnings
  are suppressed so the reported message is the informative one.

* `plot()` on a fitted design now validates `which` before touching the
  graphics state. Previously `par()` ran first, which opens the default device
  when none is active, so a mistyped `which` threw *and* left a device open,
  writing a stray `Rplots.pdf` in non-interactive sessions.

* `fit_distributions()` no longer leaves the caller's random-number stream
  altered. When a stratification variable is longer than `max_fit_n` the
  column is subsampled under a fixed seed, so that AIC family selection stays
  stable between runs; that seed is now saved and restored rather than left in
  place. Two consequences are fixed. `fit_distributions()` no longer changes
  the results of whatever the caller does next; and because
  `stratify_multi(distributions = "auto")` calls it before drawing its
  multi-start boundaries, the user's own `set.seed()` is no longer silently
  overridden for the rest of the optimisation, previously two different user
  seeds could produce byte-identical auto-mode designs, so a design looked
  reproducible while not in fact being under the user's control. Family
  selection and the fitted marginal parameters remain deterministic; pass
  `max_fit_n >= nrow(data)` to make the subsample follow the caller's seed
  instead.

* The application's *R Code* tab now emits the `set.seed()` call for the run it
  is describing, so the generated script reproduces the design shown on screen
  rather than a different draw from the multi-start scheme. The seed is
  settable in the sidebar and changing it marks the displayed results stale.

## Documentation

* `stratify_multi()` and `stratify_multi_density()` gained a
  **Reproducibility** section stating that the multi-start scheme sets no seed
  and that `set.seed()` should be used for any design that will be reported,
  with an indication of how the spread across starts shrinks as `n_starts`
  grows.

## Interactive application

* The companion 'shiny' application now ships **inside the package** under
  `inst/shiny/` and is launched with the new exported function `run_app()`.
  It exposes the full workflow (data upload/simulation/density modes, all
  weight modes and constraints, both solvers, interactive visualisations,
  benchmarking, `find_n()`, `choose_H()`, an omega explorer, a downloadable
  report, and a reproducible-R-script generator). The user-interface packages
  (`shiny`, `bslib`, `DT`, `jsonlite`, and optional `plotly`, `readxl`,
  `pagedown`) are declared in `Suggests` and used conditionally, so the core
  package still installs and checks without them.

## Reproducibility

* `REGENERATE_RESULTS.R` now also runs Scenarios C (Weibull × Pareto),
  D (trivariate, 3×3×3), and E (Gaussian copula), reports an auto-mode
  family-recovery table, and emits three publication figures as PDFs into
  `paper-figures/`: the Scenario A stratification (3-D joint density with
  boundary walls + top-view boundary grid), the parametric-vs-empirical
  convergence trajectory, and a cross-scenario relative-efficiency /
  optimality-gap summary.

## Dependencies

* `survey` added to Suggests (used by `as_svydesign()`).
