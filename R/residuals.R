#' @name residuals.gctsc
#' @title Randomized Quantile Residuals for Copula Count Time Series Models
#'
#' @description
#' Computes randomized quantile residuals for a fitted Gaussian or
#' Student--t copula count time series model.
#'
#' For discrete responses, residuals are constructed using the
#' randomized probability integral transform (PIT) as proposed by
#' Dunn and Smyth (1996). When the likelihood is evaluated via
#' simulation (TMET or GHK), the same engine is used to approximate
#' the required conditional probabilities.
#'
#' @param object A fitted model object of class \code{"gctsc"},
#'   as returned by \code{\link{gctsc}}.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @return A list of class \code{"gctsc.residuals"} containing:
#' \itemize{
#'   \item \code{residuals}: Numeric vector of randomized quantile residuals.
#'   \item \code{pit}: Numeric vector of probability integral transform values.
#' }
#'
#' @details
#' For observation \eqn{y_t}, let \eqn{F_t(y_t^-)} and \eqn{F_t(y_t)}
#' denote the conditional CDF evaluated at \eqn{y_t - 1} and \eqn{y_t},
#' respectively. The PIT value is computed as
#' \deqn{
#' e_t = F_t(y_t^-) + u_t \{F_t(y_t) - F_t(y_t^-)\},
#' }
#' where \eqn{u_t \sim \mathrm{Uniform}(0,1)}.
#'
#'For Gaussian copulas, residuals are obtained as
#'\eqn{r_t = \Phi^{-1}(e_t)}.
#'
#'For Student--t copulas with degrees of freedom \code{df},
#'the residuals are defined as \eqn{r_t = t_{\nu}^{-1}(e_t)},
#'where \eqn{t_{\nu}^{-1}} denotes the quantile function of the
#'Student--t distribution.
#' 
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996),
#' Randomized quantile residuals,
#' \emph{Journal of Computational and Graphical Statistics},
#' \strong{5}(3): 236--244.
#'
#' Nguyen, Q. N., and De Oliveira, V. (2026),
#' Likelihood Inference in Gaussian Copula Models for Count Time Series
#' via Minimax Exponential Tilting,
#' \emph{Computational Statistics and Data Analysis}.
#'
#' Nguyen, Q. N., and De Oliveira, V. (2026),
#' Scalable Likelihood Inference for Student--\eqn{t} Copula Count Time Series,
#' Manuscript in preparation.
#'
#' @examples
#' # Simulate Poisson AR(1) data under a Gaussian copula
#' set.seed(1)
#' y <- sim_poisson(mu = 5, tau = 0.7,
#'                  arma_order = c(1, 0),
#'                  nsim = 500,
#'                  family = "gaussian")$y
#'
#' fit <- gctsc(
#'   y ~ 1,
#'   data = data.frame(y = y),
#'   marginal = poisson.marg(),
#'   cormat = arma.cormat(1, 0),
#'   family = "gaussian",
#'   method = "CE",
#'   options = gctsc.opts(seed = 1, M = 1000)
#' )
#'
#' res <- residuals(fit)
#' hist(res$residuals, main = "Randomized Quantile Residuals")
#' hist(res$pit, main = "PIT Histogram")
#'
#' @importFrom stats resid
#' @export
residuals.gctsc <- function(object, ...) {
  is.int <- TRUE

  if (is.int && exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    set.seed(object$options$seed)
  }

  bounds <- object$marginal$bounds
  family <- object$family
  ab <- bounds(object$y,
               object$x,
               object$coef[object$ibeta], family = family, df= object$df)
  
  
  method <- object$method
  
  
  cfg <- list(
    method = method,
    arg2 = if (is.int) object$options$M else object$c,
    ret_llk = FALSE,
    pm = object$pm,
    od=object$cormat$od,
    QMC=object$QMC,
    df= object$df
  )


  tau <- object$coef[object$itau]
  QMC <- object$QMC
  res <- llk.fn(cfg, ab, tau,family )

  u <- runif(object$n)
  pit_vals <- rep(NA, object$n)
  res_vals <- rep(NA, object$n)


  if (family =="gaussian"){
      pit_vals <- res[,2] - u * ( res[,2] -res[,1])
      res_vals <- qnorm(pit_vals)
   
  } else {
    if (is.int) {
      pit_vals <- res[,2] - u * ( res[,2] -res[,1])
      res_vals <- qt(pit_vals,df = object$df)
    }
  }
  
  out <- list(residuals = res_vals, pit = pit_vals)
  class(out) <- "gctsc.residuals"
  return(out)
}


#' @title Diagnostic Plots for Fitted Copula Count Time Series Models
#'
#' @description
#' Produces diagnostic plots for a fitted Gaussian or Student--t copula
#' count time series model of class \code{"gctsc"}.
#'
#' The diagnostics are based on randomized quantile residuals and
#' probability integral transform (PIT) values.
#'
#' @param x A fitted model object of class \code{"gctsc"}.
#' @param caption Optional character vector of length 5 providing captions
#'   for the plots.
#' @param main Optional main titles for the plots.
#' @param level Confidence level for the Q--Q envelope (default 0.95).
#' @param col.lines Color used for reference lines.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @details
#' The following diagnostic plots are produced:
#' \enumerate{
#'   \item Time series of randomized quantile residuals.
#'   \item Q--Q plot against the reference distribution.
#'   \item Histogram of PIT values.
#'   \item Autocorrelation function (ACF) of residuals.
#'   \item Partial autocorrelation function (PACF) of residuals.
#' }
#'
#' For Gaussian copulas, residuals are compared against the standard
#' normal distribution. For Student--t copulas, residuals are compared
#' against a Student--t distribution with degrees of freedom obtained from fitted model.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @seealso \code{\link{residuals.gctsc}}
#' @export
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
#'   cormat = arma.cormat(p = 1, q = 0), family ="gaussian",
#'   method = "CE",
#'   options = gctsc.opts(seed = 1, M = 1000),
#'   c = 0.5
#' )
#'
#' # Produce diagnostic plots
#' par(mfrow = c(2, 3), ask = FALSE)
#' plot(fit)

#' @seealso \code{\link{residuals.gctsc}} for computing the residuals used in the plots.
#'
#' @export
plot.gctsc <- function(x,
                       caption = rep("", 5),
                       main = rep("", 5),
                       level = 0.95,
                       col.lines = "gray",
                       ...) {
  
  if (!inherits(x, "gctsc"))
    stop("plot.gctsc() is only for objects of class 'gctsc'.")
  
  object  <- x
  family  <- object$family
  res_out <- residuals(object)
  res     <- res_out$residuals
  pit     <- res_out$pit
  
  has_res_na <- anyNA(res)
  has_pit_na <- anyNA(pit)
  
  op <- par(mfrow = c(2, 3))
  on.exit(par(op))
  
  ## 1. Time series of residuals
  if (!has_res_na) {
    plot(res, type = "l",
         xlab = "Time",
         ylab = "Quantile residual",
         main = main[1], ...)
    abline(h = 0, col = col.lines)
    mtext(caption[1], 3, 0.25)
  } else {
    plot.new()
  }
  
  ## 2. Q-Q plot
  if (!has_res_na) {
    if (family == "t") {
      df <- object$df
      n  <- length(res)
      
      p  <- ppoints(n)
      theo <- qt(p, df = df)
      res_sorted <- sort(res)
      
      plot(theo, res_sorted,
           xlab = "Theoretical t quantiles",
           ylab = "Sorted residuals",
           main = main[2], ...)
      abline(0, 1, col = col.lines)
      
    } else {
      qqnorm(res,
             main = main[2],
             ylab = "Sorted residuals",
             ...)
      qqline(res, col = col.lines)
    }
    mtext(caption[2], 3, 0.25)
  } else {
    plot.new()
  }
  
  ## 3. PIT histogram
  if (!has_pit_na) {
    hist(pit,
         breaks = 20,
         col = "skyblue",
         border = "white",
         freq = FALSE,
         xlim = c(0, 1),
         main = main[3],
         xlab = "PIT values")
    abline(h = 1, col = col.lines, lty = 2)
    mtext(caption[3], 3, 0.25)
  } else {
    plot.new()
  }
  
  ## 4. ACF
  if (!has_res_na) {
    acf(res,
        main = main[4],
        ci.col = col.lines,
        na.action = na.pass)
    mtext(caption[4], 3, 0.25)
  } else {
    plot.new()
  }
  
  ## 5. PACF
  if (!has_res_na) {
    pacf(res,
         main = main[5],
         ci.col = col.lines,
         na.action = na.pass)
    mtext(caption[5], 3, 0.25)
  } else {
    plot.new()
  }
  
  ## 6th panel left empty
  plot.new()
  
  invisible(NULL)
}