#' Approximate Log-Likelihood via GHK Simulation
#'
#' Computes the approximate log-likelihood for a count time series model
#' based on a Gaussian or Student--t copula using the
#' Geweke--Hajivassiliou--Keane (GHK) simulator.
#'
#' The GHK method approximates the multivariate normal or Student--t
#' rectangle probability defining the copula likelihood by sequential
#' simulation from truncated conditional distributions.
#'
#' Two copula families are supported:
#' \itemize{
#'   \item \strong{Gaussian copula} via \code{pmvn_ghk()}
#'   \item \strong{Student--t copula} via \code{pmvt_ghk()}
#' }
#'
#' In both cases, the latent dependence structure is parameterized
#' through an ARMA(\eqn{p,q}) process.
#'
#' @param lower Numeric vector of length \code{n} giving the lower
#'   bounds of the transformed latent variables.
#'
#' @param upper Numeric vector of length \code{n} giving the upper
#'   bounds of the transformed latent variables.
#'
#' @param tau Numeric vector of ARMA dependence parameters ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}, where
#'   \eqn{\phi_i} are autoregressive (AR) coefficients and
#'   \eqn{\theta_j} are moving-average (MA) coefficients.
#'
#' @param od Integer vector \code{c(p, q)} specifying the AR and MA orders.
#'
#' @param M Positive integer specifying the number of Monte Carlo or
#'   quasi-Monte Carlo samples used in the simulation.
#'
#' @param QMC Logical; if \code{TRUE} (default), quasi-Monte Carlo
#'   integration is used. Otherwise, standard Monte Carlo sampling
#'   is applied.
#'
#' @param ret_llk Logical; if \code{TRUE} (default), returns the approximate
#' log-likelihood. If \code{FALSE}, internal diagnostic quantities
#' from the GHK simulator are returned. This option is primarily
#' intended for internal use and methodological research.
#'
#' @param df Degrees of freedom for the t copula. Must be greater than 2.
#'   Required only for \code{pmvt_ghk()}.
#' @param engine Character string specifying the conditional simulation
#'   engine used in the GHK approximation for the t copula.
#'   The default \code{"mvmn"} uses the multivariate mixture of normal (Standard implementation).
#'   The alternative \code{"mvt"} uses the direct Student--t
#'   conditional simulation scheme included for reproducibility of
#'   results reported in Nguyen and De Oliveira (2026).
#'   Ignored for the Gaussian copula.
#'
#' @return By default, a numeric scalar giving the approximate
#'   log-likelihood. If \code{ret_llk = FALSE}, internal diagnostic
#'   quantities from the GHK simulator are returned (primarily for
#'   internal use).
#'
#' @details
#' The GHK simulator approximates the multivariate normal or t
#' rectangle probability by decomposing it into a sequence of
#' one-dimensional truncated conditional distributions.
#'
#' @seealso \code{\link{pmvn_ghk}}, \code{\link{pmvt_ghk}}
#' @references
#' Nguyen, Q. N., & De Oliveira, V. (2026). Likelihood Inference in Gaussian Copula Models for Count Time Series
#' via Minimax Exponential Tilting
#' \emph{Journal of Computational Statistics and Data Analysis}.
#' 
#' Nguyen, Q. N., & De Oliveira, V. (2026).
#' Scalable Likelihood Inference for Student--\eqn{t} Copula Count Time Series.
#' Manuscript in preparation.
#'
#' @examples
#' ## Gaussian copula example
#' mu <- 10
#' tau <- 0.2
#' arma_order <- c(1, 0)
#'
#' sim_data <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                         nsim = 500, family = "gaussian", seed = 1)
#'
#' y <- sim_data$y
#' a <- qnorm(ppois(y - 1, lambda = mu))
#' b <- qnorm(ppois(y, lambda = mu))
#'
#' llk_gauss <- pmvn_ghk(lower = a, upper = b, tau = tau, od = arma_order,
#'                       M = 1000)
#'
#'
#' ## Student--t copula example
#' df <- 8
#'
#' sim_data_t <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                           nsim = 500, family = "t", df = df, seed = 1)
#'
#' y_t <- sim_data_t$y
#' a_t <- qt(ppois(y_t - 1, lambda = mu), df = df)
#' b_t <- qt(ppois(y_t, lambda = mu), df = df)
#'
#' llk_t <- pmvt_ghk(lower = a_t, upper = b_t, tau = tau, od = arma_order,
#'                   M = 1000, df = df)
#' @rdname pmv_ghk
#' @export
pmvn_ghk <- function(lower, upper, tau, od, M = 1000, QMC = TRUE, ret_llk = TRUE) {
  
  ghk_core(lower, upper, tau, od,
           family = "gaussian",
           M = M, QMC = QMC,
           ret_llk = ret_llk)
}

#' @rdname pmv_ghk
#' @export
pmvt_ghk <- function(lower, upper, tau, od, M = 1000, QMC = TRUE, 
                     ret_llk = TRUE, df, engine = c("mvmn", "mvt")) {
  
  engine <- match.arg(engine)
  
  ghk_core(lower, upper, tau, od,
           family = "t",
           M = M, QMC = QMC,
           ret_llk = ret_llk,
           df = df,
           engine = engine)
}


#' @keywords internal
#' @noRd
loglik_ghk <- function(ab, tau, od,
                       family,
                       M = 1000,
                       QMC = TRUE,
                       ret_llk = TRUE,
                       df = NULL,
                       engine = "mvmn") {
  
  if (length(tau) != sum(od))
    stop("Length of 'tau' must equal p+q.")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(is.na(ab)))
    return(-1e20)
  
  a <- ab[, 1]
  b <- ab[, 2]
  
  tryCatch(
    ghk_core(a, b, tau, od,
             family = family,
             M = M, QMC = QMC,
             ret_llk = ret_llk,
             df = df,
             engine = engine),
    error = function(e) {
      message("GHK failed: ", e$message)
      -1e20
    }
  )
}



ghk_core <- function(lower, upper, tau, od,
                     family,
                     M = 1000,
                     QMC = TRUE,
                     ret_llk = TRUE,
                     df = NULL,
                     engine = "mvmn") {
  
  if (any(upper < lower))
    stop("Invalid bounds: some upper bounds are smaller than lower bounds.")
  
  n <- length(lower)
  
  ghk_obj <- cond_mv_ghk(n, tau, od)
  
  # choose simulator
  if (family == "gaussian") {
    
    exp_psi <- sample_mvn(
      ghk_obj, lower, upper,
      M = M, QMC = QMC,
      ret_llk = ret_llk,
      method = "GHK"
    )
    
  } else {  # t copula
    
    method_flag <- if (engine == "mvt") "GHK_MVT" else "GHK"
    
    exp_psi <- sample_mvt(
      ghk_obj, lower, upper,
      M = M, QMC = QMC,
      ret_llk = ret_llk,
      method = method_flag,
      df = df
    )
  }
  
  if (ret_llk) {
    exponent <- exp_psi[[2]]
    return(exponent + log(exp_psi[[1]]))
  } else {
    return(exp_psi$summary_stats)
  }
}