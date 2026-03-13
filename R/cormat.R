#' ARMA Correlation Structure for Copula Count Time Series Models
#'
#' Constructs an ARMA(\eqn{p,q}) correlation structure for use in
#' Gaussian and Student--t copula count time series models.
#'
#' The ARMA model specifies the dependence structure of the latent
#' copula process. 
#'
#' @param p Non-negative integer specifying the autoregressive (AR) order.
#' @param q Non-negative integer specifying the moving-average (MA) order.
#'   The model ARMA(0,0) is not supported.
#'
#' @param tau.lower Optional numeric vector of length \code{p + q}
#'   specifying lower bounds for the ARMA parameters. 
#'
#' @param tau.upper Optional numeric vector of length \code{p + q}
#'   specifying upper bounds for the ARMA parameters.
#'
#' @return An object of class \code{"arma.gctsc"} and \code{"cormat.gctsc"}
#' containing:
#' \itemize{
#'   \item \code{npar}: Number of ARMA parameters (\eqn{p + q}).
#'   \item \code{od}: Integer vector \code{c(p, q)}.
#'   \item \code{start}: Function to compute starting values from data,
#'         typically using \code{\link[stats]{arima}}.
#'   \item \code{lower}, \code{upper}: Parameter bounds.
#' }
#'
#' @details
#' The ARMA parameters must define a stationary and invertible process.
#' These conditions are enforced during model fitting.
#'
#' @seealso \code{\link{gctsc}},
#'   \code{\link{poisson.marg}},
#'   \code{\link{predict.gctsc}}
#'
#' @export

arma.cormat <- function(p = 0, q = 0, tau.lower = NULL, tau.upper = NULL) {
  ans <- list()

  ans$npar <- p + q

  ans$start <- function(y) {
    tau <- arima(y, order = c(p, 0, q))$coef[1:(p + q)]
    ar_names <- if (p > 0) paste0("ar", 1:p) else NULL
    ma_names <- if (q > 0) paste0("ma", 1:q) else NULL
    names(tau) <- c(ar_names, ma_names)

    # Attach bounds if provided
    if (!is.null(tau.lower)) attr(tau, "lower") <- tau.lower
    if (!is.null(tau.upper)) attr(tau, "upper") <- tau.upper

    tau
  }

  ans$od <- c(p, q)
  class(ans) <- c("arma.gctsc", "cormat.gctsc")
  ans
}
