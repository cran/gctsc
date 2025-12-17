
#' Simulate latent ARMA Gaussian process
#'
#' @name sim_latent
#' @keywords internal
#' @noRd
sim_latent <- function(a, sigma2, n = 100) {
  p <- length(a$phi)
  q <- length(a$theta)
  burn1 <- 1000  #' Burn-in for AR
  burn2 <- 100   #' Burn-in for MA
  total_n <- n + burn1 + burn2
  
  if (all(a$phi == 0) && all(a$theta == 0)) {
    return(rnorm(n, mean = 0, sd = sqrt(sigma2)))
  }
  
  z <- rnorm(total_n, mean = 0, sd = sqrt(sigma2))
  
  if (q > 0) {
    z <- stats::filter(z, c(1, a$theta), sides = 1)
    z <- z[-(1:burn2)]
  }
  
  if (p > 0) {
    z <- stats::filter(z, a$phi, method = "recursive")
  }
  
  z_final <- z[(burn1 + 1):(burn1 + n)]
  stopifnot(length(z_final) == n)
  return(z_final)
}


#' Inverse CDF simulation for copula-based marginals
#'
#' @name simulate_marginal
#' @keywords internal
#' @noRd
simulate_marginal <- function(marg, u, lambda) {
  marg <- match.arg(marg, choices = c("poisson", "negbin", "zip", "binom", "zib", "bbinom", "zibb"))
  
  switch(marg,
         "poisson" = {
           mu <- lambda$mu
           qpois(u, lambda = mu)
         },
         "negbin" = {
           mu <- lambda$mu
           size <- 1 / lambda$dispersion
           qnbinom(u, mu = mu, size = size)
         },
         "zip" = {
           mu <- lambda$mu
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)    #remove the effect of u <= pi0
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qpois(v, lambda = mu))
         }
         ,
         "binom" = {
           prob <- lambda$prob
           size <- lambda$size
           qbinom(u, size, prob)
         }
         ,
         "zib" = {
           size <- lambda$size
           prob <- lambda$prob
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qbinom(v, size,prob))
         }
         ,
         "bbinom" = {
           n <- length(u)
           if (!requireNamespace("VGAM", quietly = TRUE)) {
             stop("Please install the 'VGAM' package for beta-binomial simulation.")
           }
           prob <- .recyclen(lambda$prob, n, "prob")
           rho <- lambda$rho
           size <- lambda$size
           qbbinom_custom(u, size, prob, rho)
         },
         "zibb" = {
           if (!requireNamespace("VGAM", quietly = TRUE)) {
             stop("Please install the 'VGAM' package for beta-binomial simulation.")
           }
           n <- length(u)
           size <- lambda$size
           prob <- .recyclen(lambda$prob, n, "prob")
           rho <- lambda$rho
           pi0 <- lambda$pi0
           v <- (u - pi0) / (1 - pi0)
           v <- pmax(0, pmin(1, v))
           ifelse(u <= pi0, 0, qbbinom_custom(v, size,prob, rho))
         },
         stop("Unknown marginal distribution: ", marg)
  )
}


#' Quantile function for beta-binomial using inversion
#'
#' @name qbbinom_custom
#' @keywords internal
#' @noRd
qbbinom_custom <- function(p, size, prob, rho) {
  alpha <- prob * (1 - rho) / rho
  beta <- (1 - prob) * (1 - rho) / rho
  
  sapply(seq_along(p), function(i) {
    support <- 0:size
    cdf_vals <- VGAM::pbetabinom.ab(support, size = size,
                                    shape1 = alpha[i], shape2 = beta[i])
    idx <- which(cdf_vals >= p[i])[1]
    if (is.na(idx)) size else support[idx]
  })
}

