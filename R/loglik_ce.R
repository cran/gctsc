#' Approximate Log-Likelihood via Continuous Extension (CE)
#'
#' Computes the approximate log-likelihood of a count time series model using the
#' Continuous Extension (CE) method. This approach approximates the probability
#' of the observed count vector by replacing the discrete indicator function with
#' a smooth approximation, controlled by a smoothing parameter \code{c}.
#'
#' The method is applicable to Gaussian copula models with ARMA dependence and
#' arbitrary discrete marginal distributions, provided the marginal CDF bounds
#' are given for each observation.
#'
#' @param lower Numeric vector of length \code{n}; lower bounds of the transformed latent variables.
#' @param upper Numeric vector of length \code{n}; upper bounds of the transformed latent variables.
#' @param tau Numeric vector of ARMA dependence parameters (concatenated \code{phi}, \code{theta}).
#' @param od Integer vector of length 2: \code{c(p, q)} for AR and MA orders.
#' @param c Smoothing bandwidth parameter; default is \code{0.5}. Must lie in (0,1).
#' @param ret_llk Logical; return log-likelihood if TRUE.
#'
#' @return Returns a numeric value representing the approximate log-likelihood.
#'
#' @examples
#' # Simulate Poisson AR(1) data
#' mu=10
#' tau=0.2
#' arma_order=c(1,0)
#' sim_data <- sim_poisson(mu =mu, tau=tau, arma_order=arma_order, nsim = 1000, seed = 1)
#' y <- sim_data$y
#'
#' # Compute latent bounds for CE method
#' a <- qnorm(ppois(y - 1, lambda = mu))  # lower bound
#' b <- qnorm(ppois(y, lambda = mu))      # upper bound
#'
#' # Approximate log-likelihood with CE method
#' llk_ce <- pmvn_ce(lower = a, upper = b, tau = tau, od = arma_order, c = 0.5)
#' print(llk_ce)

#'
#' @export
pmvn_ce <-  function(lower, upper, tau, od, c=0.5,  ret_llk=TRUE){
  EPS <- sqrt(.Machine$double.eps)
  EPS1 <- 1 - EPS
  if (anyNA(lower) || anyNA(upper) || any(upper < lower - EPS)) {
    bad_idx <- which(is.na(lower) | is.na(upper) | upper < lower - EPS)
    warning(
      "CE approximation failed due to invalid bounds at indices: ",
      paste(bad_idx, collapse = ", "), ". ",
      "Check for missing values or upper < lower."
    )
    return(-1e20)
  }

  n <- length(lower)
  cdf <- pnorm(upper)
  pdf <- cdf - pnorm(lower)

  # Smoothed latent quantile
  r <- qnorm(pmin(EPS1, pmax(EPS, cdf - c * pdf)))

  # Extract ARMA orders and coefficients
  .p <- od[1]; .q <- od[2]
  iar <- if (.p) 1:.p else NULL
  ima <- if (.q) (.p + 1):(.p + .q) else NULL
  phi <- tau[iar]
  theta <- tau[ima]


  # Adjust orders and coefficients if p = 0 or q = 0
  if (.p == 0) { p <- 1; phi <- 0 } else p <- .p
  if (.q == 0) { q <- 1; theta <- 0 } else q <- .q

  m <- max(p, q)
  Tau <- list(phi = phi, theta = theta)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))


  if (ret_llk) {
    model <- list(
      phi = phi, theta = theta, r = r, theta_r = theta_r,
      n = n, p = p, q = q, m = m, sigma2 = sigma2, a = pdf
    )
    results <- CE_core_recursive(gamma, model)
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


#' @keywords internal
#' @noRd
loglik_ce<- function(ab, tau, c = 0.5, od, ret_llk = TRUE) {
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported. Please specify at least one AR or MA term.")
  }

  result <- tryCatch({
    pmvn_ce(lower = ab[,1], upper=ab[,2],  tau=tau,
                   od=od, c=c,  ret_llk=ret_llk)}, error = function(e) {
                     message("CE failed: ", e$message)
                     return( -1e20)
                   })
  return(result)
}





