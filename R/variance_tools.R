# =============================================================================
# variance_tools.R
# Functions for computing and comparing composite variances, efficiency ratios,
# and CV diagnostics for stratified multivariate survey designs.
#
# Fix in this revision:
#   variance_lower_bound() computed A_h via
#   sqrt(rowSums(sweep(sigma2, 2, omega * W_h^2, "*"))), which misaligns
#   omega/W_h whenever H != q.  Now uses the corrected .compute_A_h()
#   (A_h = W_h * sqrt(sum_k omega_k sigma2_hk)) shared with neyman_alloc().
# =============================================================================


#' Compute the Composite Variance for a Stratified Design
#'
#' Evaluates \eqn{\Phi(\mathbf{B}, \mathbf{n}) = \sum_k \omega_k V_k(\mathbf{B},
#' \mathbf{n})} for a given boundary configuration and allocation vector.
#'
#' @param sigma2  Matrix of within-stratum variances (\eqn{H \times q}).
#' @param W_h     Stratum weight vector (length \eqn{H}).
#' @param n_alloc Allocated stratum sample sizes (length \eqn{H}).
#' @param N_h     Stratum population sizes (length \eqn{H}).
#' @param omega   Study-variable weight vector (length \eqn{q}).
#'
#' @return Scalar \eqn{\Phi}.
#'
#' @examples
#' sigma2 <- matrix(c(4, 9, 1, 4), nrow = 2, ncol = 2)
#' composite_variance(sigma2, W_h = c(0.4, 0.6), n_alloc = c(80, 120),
#'                   N_h = c(400, 600), omega = c(0.5, 0.5))
#'
#' @export
composite_variance <- function(sigma2, W_h, n_alloc, N_h, omega) {
  sigma2 <- as.matrix(sigma2)
  fpc    <- pmax(1 / pmax(n_alloc, 1) - 1 / pmax(N_h, 1), 0)
  sum(omega * colSums(W_h^2 * sigma2 * fpc))
}


#' Coefficients of Variation for Each Study Variable
#'
#' Returns \eqn{\mathrm{CV}(\hat{\mu}_k) = \sqrt{V_k} / \mu_k} for each study
#' variable under a given stratified design.
#'
#' @param sigma2  Matrix of within-stratum variances (\eqn{H \times q}).
#' @param mu      Matrix of within-stratum means (\eqn{H \times q}).
#' @param W_h     Stratum weight vector (length \eqn{H}).
#' @param n_alloc Stratum sample sizes (length \eqn{H}).
#' @param N_h     Stratum population sizes (length \eqn{H}).
#'
#' @return Named numeric vector of length \eqn{q}.
#'
#' @examples
#' sigma2 <- matrix(c(4, 9, 1, 4), nrow = 2, ncol = 2)
#' mu     <- matrix(c(10, 20, 5, 15), nrow = 2, ncol = 2)
#' design_cv(sigma2, mu, W_h = c(0.4, 0.6), n_alloc = c(80, 120),
#'           N_h = c(400, 600))
#'
#' @export
design_cv <- function(sigma2, mu, W_h, n_alloc, N_h) {
  sigma2 <- as.matrix(sigma2);  mu <- as.matrix(mu);  q <- ncol(sigma2)
  cvs    <- numeric(q)
  for (k in seq_len(q)) {
    fpc       <- pmax(1 / pmax(n_alloc, 1) - 1 / pmax(N_h, 1), 0)
    V_k       <- sum(W_h^2 * sigma2[, k] * fpc)
    mu_pop_k  <- sum(W_h * mu[, k])
    cvs[k]    <- if (abs(mu_pop_k) < 1e-15) NA_real_ else sqrt(V_k) / abs(mu_pop_k)
  }
  setNames(cvs, paste0("CV_Y", seq_len(q)))
}


#' Relative Efficiency of Two Stratified Designs
#'
#' Returns \eqn{\mathrm{RE} = \Phi_B / \Phi_A}; values > 1 mean design A is
#' more efficient.
#'
#' @param phi_A Composite variance of design A (scalar, must be > 0).
#' @param phi_B Composite variance of design B (scalar).
#'
#' @return Scalar relative efficiency.
#'
#' @examples
#' relative_efficiency(phi_A = 0.0042, phi_B = 0.0068)
#'
#' @export
relative_efficiency <- function(phi_A, phi_B) {
  if (phi_A <= 0) stop("'phi_A' must be positive.")
  phi_B / phi_A
}


#' Cauchy-Schwarz Variance Lower Bound
#'
#' \deqn{\Phi_{\mathrm{lb}}(\mathbf{B}) =
#'   \frac{(\sum_{\mathbf{h}} A_{\mathbf{h}})^2}{n_{\mathrm{tot}}}
#'   - \sum_{\mathbf{h},k}
#'     \frac{\omega_k W_{\mathbf{h}}^2 \sigma^2_{\mathbf{h}k}}{N_{\mathbf{h}}}}
#'
#' The bound is tight if and only if \eqn{\mathbf{n}} is the Neyman allocation.
#'
#' @param sigma2  Matrix of within-stratum variances (\eqn{H \times q}).
#' @param W_h     Stratum weight vector (length \eqn{H}).
#' @param N_h     Stratum population sizes (length \eqn{H}).
#' @param omega   Study-variable weight vector (length \eqn{q}).
#' @param n_total Total sample size.
#'
#' @return Scalar lower bound.
#'
#' @export
variance_lower_bound <- function(sigma2, W_h, N_h, omega, n_total) {
  sigma2 <- as.matrix(sigma2);  q <- ncol(sigma2)
  A_h    <- .compute_A_h(sigma2, W_h, omega)
  lb     <- sum(A_h)^2 / n_total
  for (k in seq_len(q))
    lb <- lb - omega[k] * sum(W_h^2 * sigma2[, k] / pmax(N_h, 1))
  max(lb, 0)
}
