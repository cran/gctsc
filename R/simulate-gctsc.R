#' Simulate from Gaussian Copula Time Series Models
#'
#' These functions simulate time series data from Gaussian copula models
#' with various discrete marginals and an ARMA dependence structure.
#'
#' @param mu Mean parameter(s):
#'   \itemize{
#'     \item For Poisson, ZIP, and negative binomial: \eqn{\mu > 0}.
#'     \item Scalar or vector of length \code{nsim}.
#'   }
#' @param prob Probability parameter(s) for binomial-type marginals:
#'   \eqn{0 < p < 1}, scalar or length \code{nsim}.
#' @param tau Numeric vector of ARMA coefficients, ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}. ARMA(0,0) is not supported.
#' @param arma_order Integer vector \code{c(p, q)} giving AR and MA orders.
#' @param nsim Positive integer; number of time points to simulate.
#' @param seed Optional integer to set the random seed.
#' @param dispersion Overdispersion parameter for negative binomial marginals
#'   (\eqn{\kappa > 0} in \eqn{\mathrm{Var}(Y) = \mu + \kappa \mu^2}).
#' @param pi0 Zero-inflation probability for ZIP, ZIB, and ZIBB marginals:
#'   \eqn{0 \le \pi_0 < 1}, scalar or length \code{nsim}.
#' @param rho Intra-cluster correlation for beta-binomial and ZIBB marginals
#'   (\eqn{0 < \rho < 1} in \eqn{\mathrm{Var}(Y) = n p (1-p) [1 + (n-1)\rho]}) where \eqn{n} is the number of trials.
#' @param size Number of trials for binomial-type marginals; positive integer
#'   (scalar).
#'
#' @details
#'**Marginal types:**
#'\itemize{
#'  \item {Poisson}: Counts with mean  \eqn{\mu}.
#'  \item {Negative binomial (NB)}: Overdispersed counts with mean  \eqn{\mu} and dispersion parameter \eqn{\kappa}.
#'  \item {Binomial}: Number of successes in \eqn{n} trials with success probability \eqn{p}.
#'  \item {Beta–-binomial (BB)}: Binomial with success probability \eqn{p} following a beta distribution, allowing intra-cluster correlation  \eqn{\rho}.
#'  \item {Zero--inflated Poisson (ZIP)}: Poisson with extra probability \eqn{\pi_0} of an excess zero.
#'  \item {Zero--inflated binomial (ZIB)}: Binomial with extra probability \eqn{\pi_0} of an excess zero.
#'  \item {Zero--inflated beta–binomial (ZIBB)}: Beta–binomial with extra probability \eqn{\pi_0} of an excess zero.
#'}
#'
#' **Parameterization notes:**
#' \itemize{
#'   \item Negative binomial uses \code{dispersion} (\eqn{\kappa}) to model
#'         overdispersion: larger \code{dispersion} increases variance for a fixed mean.
#'   \item Beta--binomial and ZIBB use \code{rho} as the overdispersion parameter:
#'         \eqn{\rho} is the intra-class correlation, with \eqn{\rho \to 0}
#'         giving the binomial model.
#'   \item Zero--inflated marginals include a separate \code{pi0} parameter that
#'         controls the extra probability mass at zero.
#' }
#' 
#' @return A list with components:
#' \itemize{
#'   \item \code{y}: Simulated time series data.
#'   \item \code{z}: Latent Gaussian process values.
#'   \item \code{marginal}: Marginal distribution name.
#'   \item \code{parameters}: List of parameters used.
#'   \item \code{cormat}: ARMA structure.
#' }
#' 
#' @examples
#' # Poisson example
#' sim_poisson(mu = 10, tau = c(0.2, 0.2), arma_order = c(1, 1), nsim = 100, seed = 42)
#'
#' # Negative Binomial example
#' sim_negbin(mu = 10, dispersion = 2, tau = c(0.5, 0.5), arma_order = c(1, 1))
#'
#' # Beta-Binomial example with seasonal covariates
#' n <- 100
#' xi <- numeric(n)
#' zeta <- rnorm(n)
#' for (j in 3:n) {
#'   xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
#' }
#' prob <- plogis(0.2 + 0.3 * sin(2 * pi * (1:n) / 12) +
#'              0.5 * cos(2 * pi * (1:n) / 12) + 0.3 * xi)
#' sim_zibb(prob, rho = 1/6, pi0 = 0.2, size = 24, tau = 0.5, arma_order = c(1, 0))
#'
#' @name sim_gctsc
#' @export
sim_poisson <- function(mu, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_poisson")
  .check_mu(mu, nsim, "sim_poisson")
  sim(object = list(marg = "poisson",
                    par = list(mu = mu, arma_order = arma_order),
                    tau = tau),
      nsim = nsim, seed = seed)
}


#' @rdname sim_gctsc
#' @export
sim_negbin <- function(mu, dispersion, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_negbin")
  .check_mu(mu, nsim, "sim_negbin")
  .check_dispersion(dispersion, nsim, "sim_negbin")
  sim(object = list(marg = "negbin",
                    par = list(mu = mu, dispersion = dispersion, arma_order = arma_order),
                    tau = tau),
      nsim = nsim, seed = seed)
}

#' @rdname sim_gctsc
#' @export
sim_zip <- function(mu, pi0, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_zip")
  .check_mu(mu, nsim, "sim_zip")
  .check_pi0(pi0, nsim, "sim_zip")
  sim(object = list(marg = "zip",
                    par = list(mu = mu, pi0 = pi0, arma_order = arma_order),
                    tau = tau),
      nsim = nsim, seed = seed)
}

#' @rdname sim_gctsc
#' @export
sim_binom <- function(prob, size, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_binom")
  .check_prob(prob, nsim, "sim_binom")
  .check_size(size, "sim_binom", scalar_only = TRUE)
  sim(object = list(marg = "binom",
                    par = list(prob = prob, size = size, arma_order = arma_order),
                    tau = tau),
      nsim = nsim, seed = seed)
}

#' @rdname sim_gctsc
#' @export
sim_bbinom <- function(prob, rho, size, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_bbinom")
  .check_prob(prob, nsim, "sim_bbinom")
  .check_rho(rho, nsim, "sim_bbinom")
  .check_size(size, "sim_bbinom", scalar_only = TRUE)
  sim(object = list(marg = "bbinom",
                    par = list(prob = prob, rho = rho, size = size, arma_order = arma_order),
                    tau = tau),  nsim = nsim, seed = seed)
}


#' @rdname sim_gctsc
#' @export
sim_zib <- function(prob, pi0, size, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_zib")
  .check_prob(prob, nsim, "sim_zib")
  .check_pi0(pi0, nsim, "sim_zib")
  .check_size(size, "sim_zib", scalar_only = TRUE)
  sim(object = list(marg = "zib",
                    par = list(prob = prob, pi0 = pi0, size = size, arma_order = arma_order),
                    tau = tau),nsim = nsim, seed = seed)
}

#' @rdname sim_gctsc
#' @export
sim_zibb <- function(prob, rho, pi0, size, tau, arma_order, nsim = 100, seed = NULL) {
  .check_common(nsim, tau, arma_order, seed, "sim_zibb")
  .check_prob(prob, nsim, "sim_zibb")
  .check_rho(rho, nsim, "sim_zibb")
  .check_pi0(pi0, nsim, "sim_zibb")
  .check_size(size, "sim_zibb", scalar_only = TRUE)
  sim(object = list(marg = "zibb",
                    par = list(prob = prob, rho = rho, pi0 = pi0, size = size, arma_order = arma_order),
                    tau = tau),
      nsim = nsim, seed = seed)
}
#' Simulate data from a gctsc or specification list
#'
#' @keywords internal
#' @noRd
sim <- function(object, nsim = 100, seed = NULL, ...) {
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

  # Simulate latent Gaussian process and map to marginals
  z <- sim_latent(tau_list, sigma2, n = nsim)
  u <- pnorm(z)

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

