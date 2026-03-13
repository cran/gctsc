#' Simulate data from a gctsc or specification list
#'
#' @keywords internal
#' @noRd
sim <- function(object, nsim = 100, cop , seed = NULL, ...) {
  if (missing(object)) {
    stop("Argument 'object' is required.")
  }
  if (!is.null(seed)) set.seed(seed)
  
  # Validate components
  if (!is.list(object) || is.null(object$marg) || is.null(object$par) || is.null(object$tau)) {
    stop("object must be a list with components 'marg', 'par', and 'tau'.")
  }
  
  marg <- object$marg
  par  <- object$par
  tau  <- object$tau
  
  
  if (is.null(par$arma_order)) stop("'arma_order' must be provided in object$par.")
  if (marg %in% c("poisson","negbin","zip") && is.null(par$mu)) stop("'mu' must be provided in object$par.")
  if (marg %in% c("binom","zib","bbinom","zibb") && is.null(par$prob)) stop("'prob' must be provided in object$par.")
  
  mu         <- par$mu
  prob       <- par$prob
  od         <- par$arma_order
  dispersion <- par$dispersion
  rho        <- par$rho
  size   <- par$size
  pi0   <- par$pi0
  
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported.")
  }
  
  # Convert ARMA params
  p <- od[1]
  q <- od[2]
  iar <- if (p > 0) 1:p else NULL
  ima <- if (q > 0) (p + 1):(p + q) else NULL
  phi <- if (length(iar)) tau[iar] else numeric(0)
  theta <- if (length(ima)) tau[ima] else numeric(0)
  tau_list <- list(phi = phi, theta = theta)
  sigma2 <- 1 / sum(ma.inf(tau_list)^2)
  
  # Simulate latent process
  z <- switch(cop,
              "gaussian" = gau_latent(tau_list, sigma2, nsim),
              "t"        = t_latent(tau_list, sigma2, nsim, df = par$df))
  u <- if (cop == "t") pt(z, df = par$df) else pnorm(z)
  
  # Parameters for marginals
  
  lambda <- switch(marg,
                   "poisson" = list(mu = mu),
                   "zip"     = list(mu = mu,pi0 = pi0),
                   "negbin"  = list(mu = mu, dispersion = dispersion),
                   "binom"   = list(prob = prob, size = size),
                   "zib"     = list(prob = prob, size = size, pi0 = pi0),
                   "bbinom"  = list(prob = prob, rho = rho, size = size),
                   "zibb"    = list(prob = prob, rho = rho, size = size, pi0 = pi0),
                   stop("Invalid marginal type: ", marg)
  )
  
  y <- simulate_marginal(marg, u, lambda)
  
  return(list(
    y = y,
    z = z,
    marginal = marg,
    parameters = lambda,
    cormat = list(arma_order = od, tau = tau)
  ))
}



#' Simulate latent ARMA Gaussian process
#'
#' @name sim_latent
#' @keywords internal
#' @noRd
gau_latent <- function(a, sigma2, n = 100) {
  p <- length(a$phi)
  q <- length(a$theta)
  burnin <- 1000
  z <- arima.sim(
    model = list(
      ar = if (p > 0) a$phi else NULL,
      ma = if (q > 0) a$theta else NULL
    ),
    n  = burnin + n,
    sd = sqrt(sigma2)   # sd is standard deviation
  )
  z <- as.numeric(z[(burnin + 1):(burnin + n)])
  return(z)
}



#' Simulate latent ARMA t process
#'
#' @name sim_latent
#' @keywords internal
#' @noRd
t_latent <- function(a, sigma2, n = 100, df) {
  stopifnot(df > 0)
  burnin <- 1000
  p <- length(a$phi)
  q <- length(a$theta)
  
  z <- arima.sim(
    model = list(
      ar = if (p > 0) a$phi else NULL,
      ma = if (q > 0) a$theta else NULL
    ),
    n  = burnin + n,
    sd = sqrt(sigma2)   # sd is standard deviation
  )
  z <- as.numeric(z[(burnin + 1):(burnin + n)])
  
  S <- rchisq(1, df = df)
  v <- z / sqrt(S / df)
  
  v
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

