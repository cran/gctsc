#' Approximate Log-Likelihood via Continuous Extension (CE)
#'
#' Computes the approximate log-likelihood for a count time series model
#' based on a Gaussian or Student--t copula using the Continuous Extension (CE) method.

#'
#' The CE method replaces the discrete probability mass at each observation
#' with a smooth approximation controlled by a bandwidth parameter \code{c}.
#' This yields a tractable approximation to the multivariate rectangle
#' probability defining the likelihood.
#'
#' Two copula families are supported:
#' \itemize{
#'   \item \strong{Gaussian copula} via \code{pmvn_ce()}
#'   \item \strong{t copula} via \code{pmvt_ce()}
#' }
#'
#' In both cases, the latent dependence structure is specified through
#' an ARMA(\eqn{p,q}) model.
#'
#' @param lower Numeric vector of length \code{n} giving the lower bounds
#'   of the transformed latent variables.
#' @param upper Numeric vector of length \code{n} giving the upper bounds
#'   of the transformed latent variables.
#' @param tau Numeric vector of ARMA dependence parameters ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}.
#' @param od Integer vector \code{c(p, q)} specifying AR and MA orders.
#' @param c Smoothing bandwidth parameter in \eqn{(0,1)}. Default is \code{0.5}.
#' @param ret_llk Logical; if \code{TRUE} (default), returns the approximate
#'   log-likelihood.
#' @param df Degrees of freedom for the t copula. Required only for
#'   \code{pmvt_ce()}.
#'
#' @return A numeric value giving the approximate log-likelihood.
#'
#' @details
#' The CE approximation applies to discrete marginal distributions
#' once the corresponding latent lower and upper bounds are computed.
#' The Gaussian copula version uses the standard normal cdf,
#' while the t copula version uses the Student t cdf with degrees of
#' freedom \code{df}.
#'
#' @seealso \code{\link{pmvn_ce}}, \code{\link{pmvt_ce}}
#'
#' @references
#' Nguyen, Q. N., & De Oliveira, V. (2026).
#' Approximating Gaussian copula models for count time series:
#' Connecting the distributional transform and a continuous extension.
#' \emph{Journal of Applied Statistics}.
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
#' llk_gauss <- pmvn_ce(lower = a, upper = b,
#'                      tau = tau, od = arma_order, c = 0.5)
#'
#'
#' ## t copula example
#' df <- 8
#'
#' sim_data_t <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                           nsim = 500, family = "t", df = df, seed = 1)
#'
#' y_t <- sim_data_t$y
#' a_t <- qt(ppois(y_t - 1, lambda = mu), df = df)
#' b_t <- qt(ppois(y_t, lambda = mu), df = df)
#'
#' llk_t <- pmvt_ce(lower = a_t, upper = b_t, tau = tau, od = arma_order,
#'                  c = 0.5, df = df)

#' @rdname pmv_ce
#' @export
pmvn_ce <- function(lower, upper, tau, od,
                    c = 0.5,
                    ret_llk = TRUE) {
  
  ce_core(lower, upper, tau, od,
          family = "gaussian",
          c = c,
          ret_llk = ret_llk)
}



#' @rdname pmv_ce
#' @export
pmvt_ce <- function(lower, upper, tau, od,
                    c = 0.5,
                    ret_llk = TRUE,
                    df) {
  
  ce_core(lower, upper, tau, od,
          family = "t",
          c = c,
          ret_llk = ret_llk,
          df = df)
}


#' @keywords internal
#' @noRd
loglik_ce <- function(ab, tau, od,
                      family,
                      c = 0.5,
                      ret_llk = TRUE,
                      df = NULL) {
  
  if (length(tau) != sum(od))
    stop("Length of 'tau' must equal p+q.")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  tryCatch(
    ce_core(ab[,1], ab[,2], tau, od,
            family = family,
            c = c,
            ret_llk = ret_llk,
            df = df),
    error = function(e) -1e20
  )
}


ce_core <- function(lower, upper, tau, od,
                    family,
                    c = 0.5,
                    ret_llk = TRUE,
                    df = NULL) {
  
  EPS <- sqrt(.Machine$double.eps)
  EPS1 <- 1 - EPS
  
  if (anyNA(lower) || anyNA(upper) || any(upper < lower - EPS))
    return(-1e20)
  
  n <- length(lower)
  if(ret_llk){
  # Copula CDF and density difference
  if (family == "gaussian") {
    cdf <- pnorm(upper)
    pdf <- cdf - pnorm(lower)
    r <- qnorm(pmin(EPS1, pmax(EPS, cdf - c * pdf)))
  } else {
    cdf <- pt(upper, df)
    pdf <- cdf - pt(lower, df)
    r <- t_qt_safe(cdf - c * pdf, df = df)
  }
  }
  
  # Extract ARMA
  p <- od[1]; q <- od[2]
  phi   <- if (p > 0) tau[1:p] else 0
  theta <- if (q > 0) tau[(p + 1):(p + q)] else 0
  
  if (p == 0) p <- 1
  if (q == 0) q <- 1
  
  m <- max(p, q)
  
  Tau <- list(phi = phi, theta = theta)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma  <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))
  
  if (ret_llk) {
    
    model <- list(
      phi = phi, theta = theta, r = r,
      theta_r = theta_r, n = n,
      p = p, q = q, m = m,
      sigma2 = sigma2,
      a = pdf
    )
    
    if (family == "gaussian") {
      results <- ptmvn_ce(gamma, model)
    } else {
      model$df <- df
      results <- ptmvt_ce(gamma, model)
    }
    
    return(sum(results$llk))
    
  } else {
    cond_mv <- cond_mv_ghk(n, tau, od)
    sampler_input <- list(
      a = lower,
      b = upper,
      condSd = sqrt(cond_mv$cond_var),
      M = 1000,
      phi = phi,
      q = q,
      m = m,
      Theta = cond_mv$Theta,
      QMC = TRUE
    )
    return(rtmvn_ghk(sampler_input)$summary_stats)
  }
}


t_qt_safe <- function(p, df, eps = 1e-8) {
  # Clip probabilities
  p <- pmin(1 - eps, pmax(eps, p))
  
  if (is.infinite(df) || df >= 100) {
    # Gaussian fallback for large df
    return(qnorm(p))
  }
  
  # For extreme tails, normal approx is safer
  if (any(p < 1e-6 | p > 1 - 1e-6)) {
    return(qnorm(p))  # fast, stable fallback
  }
  
  # Otherwise use qt
  return(qt(p, df))
}
