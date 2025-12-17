#' ARMA Correlation Structure for Copula Time Series
#'
#' Constructs a correlation structure object for use in Gaussian copula time series models
#' with autoregressive moving average (ARMA) dependence.
#'
#' @param p Integer. AR order (non-negative).
#' @param q Integer. MA order (non-negative).
#' @param tau.lower Optional vector of lower bounds for the ARMA parameters.
#' @param tau.upper Optional vector of upper bounds for the ARMA parameters.
#'
#' @return An object of class \code{"arma.gctsc"} and \code{"cormat.gctsc"}, containing:
#' \item{npar}{Number of ARMA parameters.}
#' \item{od}{A length-2 vector \code{c(p, q)} giving the AR and MA order.}
#' \item{start}{Function to compute starting values from data using \code{\link[stats]{arima}}.}
#'
#' @seealso \code{\link{gctsc}}, \code{\link{poisson.marg}}, \code{\link{predict.gctsc}}
#'
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
