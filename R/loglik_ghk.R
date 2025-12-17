#' GHK Log-Likelihood Approximation
#'
#' Computes the approximate log-likelihood for a count time series model using
#' the Geweke–Hajivassiliou–Keane (GHK) simulator. This method evaluates
#' the multivariate normal rectangle probability by sequentially sampling
#' from truncated conditionals implied by the ARMA Gaussian copula.
#'
#' @param lower A numeric vector of lower truncation bounds.
#' @param upper A numeric vector of upper truncation bounds.
#' @param tau A numeric vector of ARMA dependence parameters.
#' @param od Integer vector \code{c(p, q)} specifying the AR and MA orders.
#' @param M Integer. Number of Monte Carlo or QMC samples.
#' @param QMC Logical. If \code{TRUE}, use quasi-Monte Carlo integration.
#' @param ret_llk Logical. Default is \code{TRUE} to return log-likelihood; otherwise, return diagnostic output.
#'
#' @return If \code{ret_llk = TRUE}, a numeric scalar (log-likelihood); else a list containing diagnostic statistics.
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
#' # Approximate log-likelihood with CE method
#' llk_tmet <- pmvn_ghk(lower = a, upper = b, tau = 0.2, od = c(1,0))
#' print(llk_tmet)
#'
#' @export
pmvn_ghk<- function(lower, upper, tau, od ,M = 1000, QMC = TRUE,  ret_llk=TRUE){


  n<-length(lower)
  p<- od[1]
  q <- od[2]

  if (any(upper < lower)) {
    stop("Invalid MVN bounds: some upper bounds are smaller than lower bounds.")

  }

  # compute conditional variance
  ghk_obj <- cond_mv_ghk(n,tau,od)

  # compute MVN probs and est error ---------------------------------
  exp_psi <- sample_latent(ghk_obj, lower, upper, M = M, QMC=QMC, ret_llk=ret_llk, method = "GHK")


  if(ret_llk){
    exponent <- exp_psi[[2]]
    log_pmvn <- exponent +log(exp_psi[[1]])
    return (log_pmvn)
    # weights <- exp_psi[[3]]
  }else{
    return(exp_psi$summary_stats)
  }

}

#' @keywords internal
#' @noRd
loglik_ghk <- function(ab, tau, M = 1000, od, QMC=TRUE, ret_llk = TRUE) {
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported. Please specify at least one AR or MA term.")
  }

  if (any(is.na(ab)) || any(is.nan(ab))) return(-1e20)
  a <- ab[,1]
  b <- ab[,2]
  p <- od[1]
  q <- od[2]

  result <- tryCatch({
    pmvn_ghk(a, b, tau, od,M,
                    QMC = QMC, ret_llk = ret_llk)
  }, error = function(e) {
    message("GHK failed: ", e$message)
    return(-1e20)
  })
  return(result)


}

#' @keywords internal
#' @noRd
cond_mv_ghk <- function(n, tau, od) {
  .p <- od[1]
  .q <- od[2]
  iar <- if ( .p ) 1:.p else NULL
  ima <- if ( .q ) (.p+1):(.p+.q) else NULL
  phi <- tau[iar]
  theta<-tau[ima]

  # Use default (p = 1, q = 1) if either AR or MA order is zero
  if(.p==0){
    p<- 1
    phi <-0
  } else {p<-.p}
  if(.q==0){
    q<-1
    theta<-0
  } else{q <- .q}

  m = max(p,q)
  Tau <- list(phi=phi, theta=theta)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)

  gamma <- aacvf(Tau, n - 1)

  theta_r <- c(1, theta, numeric(n))

  model <- list(phi = phi, theta_r = theta_r, p = p, q = q, m = m,
                sigma2 = sigma2, gamma = gamma, n = n)
  out_cpp <- compute_cond_var(gamma, model)
  v <- out_cpp$v
  cond_var <- v*sigma2

  ##### theta_tk is saved to compute the conditional mean
  Theta <- out_cpp$Theta

  list(
    cond_var = cond_var,
    Theta = Theta,
    m =m,
    phi = phi,
    p = .p,
    q = .q
  )
}

