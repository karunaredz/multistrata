<!-- README.md is written by hand; edit this file directly. -->

# multistrata

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/multistrata)](https://CRAN.R-project.org/package=multistrata)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->
<!-- Once a public repository exists, add its R-CMD-check badge here and
     restore URL / BugReports in DESCRIPTION. The workflow is already in
     .github/workflows/R-CMD-check.yaml. -->


**Multivariate optimum stratification and sample allocation, solved jointly.**

Stratified survey design asks two questions at once: where should the strata
boundaries fall, and how should the sample be spread across them?
`multistrata` answers both in a single optimisation, across an arbitrary
number of stratification variables and study variables, subject to
sample-size, cost and precision constraints.

Solving the two together broadens the search. A joint optimiser considers the
output of boundary-first procedures as well as designs those procedures cannot
reach. It therefore cannot perform worse under the same objective, although the
two approaches can tie. On the synthetic `pacific_agri` population shipped with
the package, the sequential benchmark has **81% more** composite variance than
the returned joint design under the settings below.

## Installation

```r
install.packages("multistrata")
```

Or from a local source tarball:

```r
install.packages("multistrata_1.0-1.tar.gz", repos = NULL, type = "source")
```

## A design in ten lines

```r
library(multistrata)
data(pacific_agri)

res <- stratify_multi(
  data       = pacific_agri[, c("land_area", "labour_units")],
  study_vars = pacific_agri[, c("crop_yield", "livestock_value", "expenditure")],
  H_vec      = c(4, 3),            # 4 x 3 = 12 strata
  n_total    = 600,
  omega      = c(0.4, 0.4, 0.2),   # priority weights on the study variables
  n_starts   = 25)

res
```

```
================================================
  Multivariate Optimum Stratification (multistrata)
================================================
  Stratification variables : 2  (land_area, labour_units)
  Study variables          : 3  (crop_yield, livestock_value, expenditure)
  Strata per variable      : 4 x 3
  Total strata             : 12
  Population size (N)      : 8000
  Total sample (n)         : 600
  Constraint type          : n_total
  Weight mode              : empirical
  Boundary solver          : COBYLA
------------------------------------------------
  Composite variance (Phi) : 0.0060
  Optimality gap           : 0.00%
  KKT residual (scaled)    : 5.40e-02  [>= tol_kkt]
  Converged                : YES
  Multi-starts used        : 25
------------------------------------------------
  Optimal strata boundaries:
    land_area: [0.1164 | 2.581 | 5.9389 | 12.1828 | 74.5662]
    labour_units: [0.0500 | 2.2553 | 3.8982 | 13.0665]
```

(Produced under `set.seed(1)`; see [Reproducibility](#reproducibility).)

How much does solving jointly buy? `benchmark_designs()` puts the joint design
against the alternatives on the same population:

```r
benchmark_designs(
  pacific_agri[, c("land_area", "labour_units")],
  pacific_agri[, c("crop_yield", "livestock_value", "expenditure")],
  H_vec = c(4, 3), n_total = 600, omega = c(0.4, 0.4, 0.2),
  result_joint = res)
```

| Design | Composite variance | Relative efficiency | Excess over joint |
|---|---|---|---|
| Joint (`multistrata`) | 0.006035 | 1.000 | n/a |
| Equal allocation | 0.008586 | 1.423 | 42.3% |
| Proportional allocation | 0.009533 | 1.580 | 58.0% |
| Sequential (quantile + Neyman) | 0.010927 | 1.811 | 81.1% |

## The interactive application

The whole workflow is available without writing code:

```r
run_app()
```

Eleven tabs cover design, visualisation, benchmarking, precision, sample-size
search, choice of strata count, the priority-weight trade-off, deployment, a
downloadable report, and a generated R script that reproduces whatever you
built. The interface packages (`shiny`, `bslib`, `DT`, `jsonlite`) are in
`Suggests`; see `vignette("multistrata-intro")` and
`system.file("shiny", "README.md", package = "multistrata")`.

## What the package does

**Design**

| Function | Purpose |
|---|---|
| `stratify_multi()` | Joint boundaries + allocation from data (empirical or parametric weights) |
| `stratify_multi_density()` | Data-free design from marginal densities alone |
| `neyman_alloc()` | Exact generalised Neyman (water-filling) allocation |
| `fit_distributions()` | AIC family selection per stratification variable |

**Design questions**

| Function | Question it answers |
|---|---|
| `find_n()` / `find_n_density()` | What is the smallest *n* meeting my CV targets? |
| `choose_H()` / `choose_H_density()` | How many strata are worth having? |
| `omega_sweep()` / `omega_sweep_density()` | What does protecting one variable cost the others? |
| `benchmark_designs()` / `benchmark_density_designs()` | How much does joint optimisation buy? |

**Diagnostics and deployment**

| Function | Purpose |
|---|---|
| `design_cv()`, `composite_variance()`, `relative_efficiency()` | Precision of a design |
| `variance_lower_bound()` | Cauchy-Schwarz allocation bound for the returned boundaries |
| `plot_strata()`, `plot_strata_rgl()` | Boundary diagrams, 1-D through *n*-D, plus interactive 3-D |
| `assign_strata()` | Score any frame into the fitted strata |
| `as_svydesign()` | Export to **survey** / **srvyr** with design weights attached |

**Features.** Empirical or parametric stratum weights; thirteen parametric
marginal families (`norm`, `lnorm`, `gamma`, `exp`, `weibull`, `beta`, `unif`,
`chisq`, `t`, `cauchy`, `logis`, `pareto`, `triang`) with automatic AIC
selection; derivative-free COBYLA or gradient-based SLSQP boundary updates;
sample-size, cost (with fixed overhead) and coefficient-of-variation
constraints; multi-start global search with a reported optimality gap and KKT
residuals.

## Example populations

| Dataset | Size | Variables |
|---|---|---|
| `pacific_agri` | 8,000 × 5 | `land_area`, `labour_units`, `crop_yield`, `livestock_value`, `expenditure` |
| `health_survey` | 6,000 × 4 | `age`, `bmi`, `systolic_bp`, `cholesterol` |
| `business_survey` | 5,000 × 4 | `employees`, `assets`, `turnover`, `wage_bill` |

All three are synthetic but realistic, generated by
`data-raw/generate_datasets.R` under a fixed seed so every example in the
package and the paper is reproducible.

## Reproducibility

`stratify_multi()` uses a multi-start scheme with randomly drawn initial
boundaries, so **set a seed** when you need a reproducible design. The spread
across starts shrinks quickly as `n_starts` grows, but it is not zero.

## Related packages

- [**stratifyR**](https://CRAN.R-project.org/package=stratifyR), optimum
  stratification for a *single* stratification variable.
- [**FitVerse**](https://CRAN.R-project.org/package=FitVerse), parametric
  distribution fitting; used as an optional alternative engine for
  `distributions = "auto"`.

## Citation

```r
citation("multistrata")
```

Reddy, K. G. and Khan, M. G. M. (2026). Simultaneous Multivariate Optimum
Stratification and Sample Allocation under Survey Design Constraints.
*Journal of Survey Statistics and Methodology*. Forthcoming.

## License

GPL-3
