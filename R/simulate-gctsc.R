#' Simulate from Gaussian and t Copula Time Series Models
#'
#' These functions simulate time series data from Gaussian and t copula models
#' with various discrete marginals and an ARMA dependence structure.
#'
#' @param mu Mean parameter(s) for Poisson-, ZIP-, and negative
#'   binomial-type marginals. Must satisfy \eqn{\mu > 0}. May be specified
#'   as a scalar or as a numeric vector of length \code{nsim} to allow
#'   time-varying means.
#'
#' @param prob Success probability parameter(s) for binomial-type marginals.
#'   Must satisfy \eqn{0 < p < 1}. May be a scalar or a numeric vector of
#'   length \code{nsim}.
#'
#' @param tau Numeric vector of ARMA dependence coefficients, ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}, where
#'   \eqn{\phi_i} are autoregressive (AR) coefficients and
#'   \eqn{\theta_j} are moving-average (MA) coefficients.
#'   The model \code{ARMA(0, 0)} is not supported.
#'
#' @param arma_order Integer vector \code{c(p, q)} specifying the AR and MA orders.
#'
#' @param nsim Positive integer giving the number of time points to simulate.
#'
#' @param seed Optional integer used to set the random seed.
#'
#' @param dispersion Overdispersion parameter for negative binomial marginals.
#'   Must satisfy \eqn{\kappa > 0}, where
#'   \eqn{\mathrm{Var}(Y) = \mu + \kappa \mu^2}.
#'   May be a scalar or a numeric vector of length \code{nsim}.
#'
#' @param pi0 Zero-inflation probability for ZIP, ZIB, and ZIBB marginals.
#'   Must satisfy \eqn{0 \le \pi_0 < 1}. May be a scalar or a numeric vector
#'   of length \code{nsim}.
#'
#' @param rho Intra-class correlation parameter for beta-binomial and ZIBB
#'   marginals. Must satisfy \eqn{0 < \rho < 1}, where
#'   \eqn{\mathrm{Var}(Y) = n p (1-p)\{1 + (n-1)\rho\}} and \eqn{n} is the
#'   number of trials. May be a scalar or a numeric vector of length
#'   \code{nsim}.
#'
#' @param size Number of trials for binomial-type marginals; a positive
#'   integer scalar.
#'
#' @param df Degrees of freedom for the t copula. Must be a single numeric
#'   value greater than 2. Required only when \code{family = "t"}.
#'
#' @param family Character string specifying the copula family:
#'   \code{"gaussian"} or \code{"t"}.
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
#' sim_poisson(mu = 10, tau = c(0.2, 0.2), 
#'   arma_order = c(1, 1), nsim = 100, 
#'   family = "gaussian", seed = 42)
#'
#' # Negative Binomial example
#' sim_negbin(mu = 10, dispersion = 2, tau = c(0.5, 0.5),
#'   arma_order = c(1, 1),family = "gaussian",
#'   nsim = 100, seed =1)
#'
#' # Zero Inflated Beta-Binomial example with seasonal covariates
#' n <- 100
#' xi <- numeric(n)
#' zeta <- rnorm(n)
#' for (j in 3:n) {
#'   xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
#' }
#' prob <- plogis(0.2 + 0.3 * sin(2 * pi * (1:n) / 12) +
#'              0.5 * cos(2 * pi * (1:n) / 12) + 0.3 * xi)
#' sim_zibb(prob, rho = 1/6, pi0 = 0.2, size = 24, tau = 0.5,
#'  arma_order = c(1, 0),family = "t", df = 10, nsim = 100)
#'
#' @name sim_gctsc
#' @export
sim_poisson <- function(mu, tau, arma_order, nsim,
                        family = c("gaussian","t"),
                        df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df, "sim_poisson")
  .check_mu(mu, nsim, "sim_poisson")
  
  par <- list(mu = mu, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "poisson", par = par, tau = tau),
      nsim = nsim, seed = seed, cop = family)
}




#' @rdname sim_gctsc
#' @export
sim_negbin <- function(mu, dispersion, tau, arma_order, nsim = 100, 
                       family = c("gaussian","t"),
                       df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_negbin")
  .check_mu(mu, nsim, "sim_negbin")
  .check_dispersion(dispersion, nsim, "sim_negbin")
  
  par = list(mu = mu, dispersion = dispersion, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "negbin",par=par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
  
}

#' @rdname sim_gctsc
#' @export
sim_zip <- function(mu, pi0, tau, arma_order, nsim = 100, family = c("gaussian","t"), 
                    df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_zip")
  .check_mu(mu, nsim, "sim_zip")
  .check_pi0(pi0, nsim, "sim_zip")
  
  par = list(mu = mu, pi0 = pi0, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zip",par = par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_binom <- function(prob, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                      df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_binom")
  .check_prob(prob, nsim, "sim_binom")
  .check_size(size, "sim_binom", scalar_only = TRUE)
  
  par = list(prob = prob, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "binom",par=par, tau = tau),
      nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_bbinom <- function(prob, rho, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                  df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_bbinom")
  .check_prob(prob, nsim, "sim_bbinom")
  .check_rho(rho, nsim, "sim_bbinom")
  .check_size(size, "sim_bbinom", scalar_only = TRUE)
  
  par = list(prob = prob, rho = rho, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "bbinom", par= par ,
                    tau = tau),  nsim = nsim, seed = seed,  cop = family)
}


#' @rdname sim_gctsc
#' @export
sim_zib <- function(prob, pi0, size, tau, arma_order, nsim = 100, family = c("gaussian","t"),
                    df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed, family , df,  "sim_zib")
  .check_prob(prob, nsim, "sim_zib")
  .check_pi0(pi0, nsim, "sim_zib")
  .check_size(size, "sim_zib", scalar_only = TRUE)
  
  par = list(prob = prob, pi0 = pi0, size = size, arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zib",
                    par = par,
                    tau = tau),nsim = nsim, seed = seed,  cop = family)
}

#' @rdname sim_gctsc
#' @export
sim_zibb <- function(prob, rho, pi0, size, tau, arma_order, nsim = 100, 
                     family = c("gaussian", "t"),  df = NULL, seed = NULL) {
  family <- match.arg(family)
  .check_common(nsim, tau, arma_order, seed,  family , df, "sim_zibb")
  .check_prob(prob, nsim, "sim_zibb")
  .check_rho(rho, nsim, "sim_zibb")
  .check_pi0(pi0, nsim, "sim_zibb")
  .check_size(size, "sim_zibb", scalar_only = TRUE)
  
  par = list(prob = prob, rho = rho, pi0 = pi0, size = size, 
             arma_order = arma_order)
  if (family == "t") par$df <- df
  
  sim(object = list(marg = "zibb",
                    par = par,
                    tau = tau), nsim = nsim, seed = seed,  cop = family)
}


