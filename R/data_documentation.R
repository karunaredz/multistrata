# =============================================================================
# data_documentation.R
#
# Documentation for the three fictitious populations shipped with the
# package.  The datasets are generated (seeded, fully reproducible) by
# data-raw/generate_datasets.R; pacific_agri is the population analysed in
# Section 6 (Real-Data Application) of Reddy & Khan (forthcoming).
# =============================================================================


#' Pacific Agricultural Survey Population (Fictitious)
#'
#' A synthetic household-level agricultural census frame mirroring the
#' structure of Pacific-region agricultural surveys discussed in
#' Reddy & Khan (2020).  The stratification variables have lognormal and
#' gamma marginals by construction, so the parametric mode of
#' \code{\link{stratify_multi}} is correctly specified for this population.
#' Used in the real-data application of Reddy & Khan (forthcoming).
#'
#' @format A \code{data.frame} with 8000 rows and 5 variables:
#' \describe{
#'   \item{land_area}{Cultivated land area (hectares); lognormal marginal.
#'     Stratification variable.}
#'   \item{labour_units}{Full-time-equivalent farm labour; gamma marginal.
#'     Stratification variable.}
#'   \item{crop_yield}{Annual crop output (tonnes).  Study variable.}
#'   \item{livestock_value}{Livestock holdings value ('000 FJD).  Study
#'     variable.}
#'   \item{expenditure}{Annual farm input expenditure ('000 FJD).  Study
#'     variable.}
#' }
#' @usage data(pacific_agri)
#' @source Fictitious data generated for illustration purposes by
#'   \code{data-raw/generate_datasets.R} (seed 20260612).
"pacific_agri"


#' Business Turnover Survey Population (Fictitious)
#'
#' A synthetic enterprise frame with a lognormal employment distribution and
#' a Pareto-tailed asset distribution, typical of establishment surveys.
#'
#' @format A \code{data.frame} with 5000 rows and 4 variables:
#' \describe{
#'   \item{employees}{Employment head count; lognormal marginal.
#'     Stratification variable.}
#'   \item{assets}{Total assets ('000); Pareto(2.8, 1.5) marginal, top 0.5\%
#'     trimmed.  Stratification variable.}
#'   \item{turnover}{Annual turnover ('000).  Study variable.}
#'   \item{wage_bill}{Annual wage bill ('000).  Study variable.}
#' }
#' @usage data(business_survey)
#' @source Fictitious data generated for illustration purposes by
#'   \code{data-raw/generate_datasets.R} (seed 20260612).
"business_survey"


#' Health Survey Population (Fictitious)
#'
#' A synthetic individual-level health frame with a beta-shaped age
#' distribution (rescaled to 18--90 years) and a normal BMI distribution.
#'
#' @format A \code{data.frame} with 6000 rows and 4 variables:
#' \describe{
#'   \item{age}{Age in years, 18--90; beta(1.8, 2.6) shape.  Stratification
#'     variable.}
#'   \item{bmi}{Body-mass index (kg/m\eqn{^2}); normal marginal.
#'     Stratification variable.}
#'   \item{systolic_bp}{Systolic blood pressure (mmHg).  Study variable.}
#'   \item{cholesterol}{Total cholesterol (mmol/L).  Study variable.}
#' }
#' @usage data(health_survey)
#' @source Fictitious data generated for illustration purposes by
#'   \code{data-raw/generate_datasets.R} (seed 20260612).
"health_survey"
