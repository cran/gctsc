#' @name residuals.gctsc
#' @title Compute Randomized Quantile Residuals for Gaussian Copula Time Series
#'
#' @description Computes residuals for a fitted \code{gctsc} object using randomized quantile residuals.
#'
#'
#' @param object A fitted object of class \code{gctsc}, produced by \code{\link{gctsc}}.
#' @param method Can be \code{TMET} or \code{GHK}
#' @param ... Ignored. Included for S3 method compatibility.
#' @return A list containing:
#'   \item{residuals}{A numeric vector of randomized quantile residuals.}
#'   \item{pit}{A numeric vector of PIT values.}
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996), Randomized quantile residuals,
#' \emph{Journal of Computational and Graphical Statistics}, \strong{5}(3): 236-244.
#'
#' @examples
#' y <- sim_poisson(mu = 5, tau = 0.7, arma_order = c(1, 0), nsim = 100)
#' fit <- gctsc(y ~ 1, marginal = poisson.marg(), cormat = arma.cormat(1, 0),
#'              method = "GHK", options = gctsc.opts(seed = 1, M = 1000))
#' res <- residuals(fit)
#' hist(res$residuals, main = "Randomized Quantile Residuals", xlab = "Residual")
#' hist(res$pit, main = "PIT Histogram", xlab = "PIT values")
#'
#' @importFrom stats resid
#' @export
residuals.gctsc <- function(object,method = NULL, ...) {
  # is.int <- (object$method != "CE")
  is.int <- TRUE

  if (!is.null(object$options$seed)) {
    set.seed(object$options$seed)
  }

  bounds <- object$marginal$bounds
  ab <- bounds(object$y,
               object$x,
               object$coef[object$ibeta])
  
  if (is.null(method)){
    method =object$method
  }
  
  cfg <- list(
    method = method,
    arg2 = if (is.int) object$options$M else object$c,
    ret_llk = FALSE,
    pm = object$pm,
    od=object$cormat$od,
    QMC=object$QMC
  )


  tau <- object$coef[object$itau]
  QMC <- object$QMC
  res <- llk.fn(cfg, ab, tau)

  u <- runif(object$n)
  pit_vals <- rep(NA, object$n)
  res_vals <- rep(NA, object$n)

  if (is.int) {
    pit_vals <- res[,2] - u * ( res[,2] -res[,1])
    res_vals <- qnorm(pit_vals)
  } else {
    pit_vals <- pnorm(res[,2])
    res_vals <- res[,2]
  }

  out <- list(residuals = res_vals, pit = pit_vals)
  class(out) <- "gctsc.residuals"
  return(out)
}


#' Diagnostic Plots for Fitted Gaussian Copula Time Series Models
#'
#' Produces a set of diagnostic plots based on randomized quantile residuals
#' and PIT values for objects of class \code{gctsc}.
#'
#' The function displays up to five plots: time series of residuals, Q--Q plot,
#' PIT histogram, and ACF/PACF plots of the residuals. These plots help assess
#' model fit and potential misspecification.
#'
#' @param x An object of class \code{gctsc}, the result of a call to \code{\link{gctsc}}.
#' @param caption Optional character vector of length 6 to use as captions for each plot.
#' @param main Optional main title for the plots (recycled if shorter than the number of plots shown).
#' @param level Confidence level for the Q--Q plot envelope (default is 0.95).
#' @param col.lines Color for reference lines in residual and ACF/PACF plots.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @details
#' The five diagnostic plots shown are:
#'
#'
#' \enumerate{
#'   \item Time series plot of randomized quantile residuals.
#'   \item Q--Q plot comparing residuals to a standard normal distribution.
#'   \item Histogram of probability integral transform (PIT) values.
#'   \item Autocorrelation function (ACF) of the residuals.
#'   \item Partial autocorrelation function (PACF) of the residuals.
#' }
#'
#' @return This function is called for its side effects and returns \code{invisible()}.
#' @examples
#' # Simulate data from a Poisson AR(1) model
#' set.seed(123)
#' n <- 2000
#' mu <- 5
#' phi <- 0.5
#' arma_order <- c(1, 0)
#' y <- sim_poisson(mu = mu, tau = phi, arma_order = arma_order, nsim = n)$y
#'
#' # Fit the model using the CE method
#' fit <- gctsc(y~1,
#'   marginal = poisson.marg(link = "identity", lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0),
#'   method = "CE",
#'   options = gctsc.opts(seed = 1, M = 1000),
#'   c = 0.5
#' )
#'
#' # Produce diagnostic plots
#' plot(fit)

#' @seealso \code{\link{residuals.gctsc}} for computing the residuals used in the plots.
#'
#' @export
plot.gctsc <- function(x, caption = rep("", 5),
                       main = rep("", 5),
                       level = 0.95, col.lines = "gray", ...) {
  object <- x
  if (!inherits(object, "gctsc"))
    stop("use only with \"gctsc\" objects")
  
  res_output <- residuals(object)
  res <- res_output$residuals
  pit <- res_output$pit
  which <- 1:5
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  Main <- rep("", 5)
  Main[which] <- rep(main, length.out = sum(show))
  
  # Flags for NA safety
  has_res_na <- anyNA(res)
  has_pit_na <- anyNA(pit)
  
  if (has_res_na) {
    message("Some residual values are NA - plots depending on residuals will be skipped.")
  }
  if (has_pit_na) {
    message("Some PIT values are NA - plots depending on PIT will be skipped.")
  }
  
  if (show[1] && !has_res_na) {
    plot(seq_along(res), res, xlab = "Obs. number", ylab = "Quantile residual",
         main = Main[1], type = "l", ...)
    abline(h = 0, col = col.lines)
    mtext(caption[1], 3, 0.25)
  }
  if (show[2] && !has_res_na) {
    car::qqPlot(res, envelope = level, grid = FALSE,
                xlab = "Normal quantiles", ylab = "Sorted quantile residuals",
                main = Main[2], col.lines = col.lines)
    mtext(caption[2], 3, 0.25)
  }
  if (show[3] && !has_pit_na) {
    hist(pit, breaks = 20, col = "skyblue", border = "white",
         main = "PIT Histogram", xlab = "PIT values", xlim = c(0, 1), freq = FALSE)
    abline(h = 1, col = "red", lty = 2)
  }
  if (show[4] && !has_res_na) {
    plot(acf(res, na.action = na.pass, plot = FALSE), ci.col = col.lines,
         main = Main[3], ...)
    mtext(caption[3], 3, 0.25)
  }
  if (show[5] && !has_res_na) {
    plot(pacf(res, na.action = na.pass, plot = FALSE), ci.col = col.lines,
         main = Main[3], ...)
    mtext(caption[3], 3, 0.25)
  }
  
  invisible()
}


