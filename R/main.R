#' Fit a Gaussian Copula Time Series Model for Count Data
#'
#' Fits a Gaussian copula model to univariate count time series using discrete
#' marginals (Poisson, negative binomial, binomial/beta–binomial, and
#' zero–inflated variants) and latent dependence through ARMA correlation
#' structures. The multivariate normal rectangle probability is evaluated using
#' Minimax Exponential Tilting (TMET), the
#' Geweke–Hajivassiliou–Keane (GHK) simulator, or the Continuous Extension (CE)
#' approximation.
#'
#' The interface mirrors \code{glm()}. Zero–inflated marginals accept a list of
#' formulas, for example \code{list(mu = y ~ x, pi0 = ~ z)}. Non–zero–inflated
#' marginals accept a single formula (e.g., \code{y ~ x1 + x2}).
#'
#' @param formula A formula (e.g., \code{y ~ x1 + x2}) or, for zero–inflated
#'   marginals, a named list of formulas \code{list(mu = ..., pi0 = ...)}.
#' @param data A data frame containing \code{y} and covariates referenced in the
#'   formula(s).
#' @param marginal A marginal model object such as \code{\link{poisson.marg}},
#'   \code{\link{negbin.marg}}, \code{\link{zib.marg}}, or \code{\link{zibb.marg}}.
#'   See \code{\link{marginal.gctsc}} for the full list of supported marginals.
#' @param cormat A correlation structure such as \code{\link{arma.cormat}}.
#' @param method One of \code{"TMET"}, \code{"GHK"}, or \code{"CE"}.
#' @param c Smoothing constant for the CE approximation (ignored for TMET/GHK).
#'   Default is \code{0.5}.
#' @param QMC Logical; use quasi–Monte Carlo sampling for simulation–based
#'   likelihood approximations.
#' @param pm Integer; truncated AR order used to approximate ARMA(\eqn{p,q})
#'   when \eqn{q>0} (TMET only).
#' @param start Optional numeric vector of starting values (marginal parameters
#'   followed by dependence parameters).
#'
#' @param options Optional list specifying the random seed and Monte Carlo
#'   sample size for the likelihood approximation. If omitted, default values
#'   from \code{gctsc.opts()} are used (i.e., \code{seed = NULL} and
#'   \code{M = 1000}).
#'
#'   For simulation–based likelihoods (\code{method = "GHK"} or
#'   \code{method = "TMET"}), providing a fixed \code{seed} (e.g.,
#'   \code{options = gctsc.opts(seed = 123)}) is **strongly recommended** to
#'   ensure reproducible and numerically stable likelihood evaluations. Without
#'   a fixed seed, the approximate log-likelihood may vary slightly between
#'   evaluations, which can affect numerical optimization.
#'
#'   See \code{\link{gctsc.opts}} for additional tuning parameters and their
#'   default values.
#'
#' @details
#'
#' **Formula handling**  
#' For zero–inflated marginals:
#' \itemize{
#'   \item If neither \code{mu} nor \code{pi0} is supplied, both default to
#'     intercept-only models (\code{mu ~ 1}, \code{pi0 ~ 1}).
#'   \item If \code{mu} is supplied but \code{pi0} is omitted,
#'     \code{pi0 ~ 1} is used.
#' }
#'
#' **Dependence structure**  
#' ARMA parameters are determined by \code{cormat}.  
#' ARMA(0,0) is not supported.
#'
#' **Method-specific notes**  
#' \itemize{
#'   \item CE ignores \code{QMC} and \code{options$M}.
#'   \item TMET and GHK require \code{options$M} to be a positive integer.
#'   \item TMET optionally uses \code{pm} when \eqn{q>0}.
#' }
#'
#' @return
#' An object of class \code{"gctsc"} containing:
#' \item{coef}{Numeric vector containing the parameter estimates, with marginal parameters first and dependence parameters last.}
#' \item{maximum}{The maximized approximate log-likelihood value, returned on the internal minus scale used by the optimizer.}
#' \item{hessian}{Estimated Hessian matrix at the optimum, when available.}
#' \item{se}{Standard errors from the inverse Hessian, or NA if unavailable or not positive-definite.}
#' \item{marginal}{The marginal model object used (e.g., Poisson, negative binomial, ZIB, ZIBB).}
#' \item{cormat}{The correlation structure object used, typically created by \code{arma.cormat()}.}
#' \item{ibeta}{Integer vector giving the indices of the marginal parameters within \code{coef}.}
#' \item{itau}{Integer vector giving the indices of the dependence parameters within \code{coef}.}
#' \item{nbeta}{Number of marginal parameters.}
#' \item{ntau}{Number of dependence parameters.}
#' \item{options}{List of fitting and optimization controls, typically created by \code{gctsc.opts()}.}
#' \item{call}{The matched call.}
#' \item{formula}{The model formula(s) used to specify the marginal mean and zero-inflation components.}
#' \item{terms}{Terms objects for the marginal model(s).}
#' \item{model}{The model frame used for estimation (returned only when needed).}
#' \item{x}{Design matrix for the marginal mean component.}
#' \item{z}{Design matrix for the zero-inflation component (ZIB or ZIBB only), or NULL otherwise.}
#' \item{y}{The response vector.}
#' \item{n}{Number of observations used in the fit.}
#' \item{method}{Character string identifying the likelihood approximation used (TMET, GHK, or CE).}
#' \item{QMC}{Logical flag indicating whether quasi-Monte Carlo sampling was used.}
#' \item{pm}{Truncated AR order used to approximate ARMA(p,q) when q > 0 (TMET only).}
#' \item{convergence}{Optimizer convergence code (0 indicates successful convergence).}
#'
#' The returned object can be used with
#' \code{\link{summary.gctsc}}, \code{\link{predict.gctsc}},
#' \code{\link{residuals.gctsc}}, and \code{\link{plot.gctsc}}.
#' 
#' @examples
#' set.seed(42)
#' n <- 200
#' y <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1,0), nsim = n)$y
#' fit <- gctsc(y ~ 1,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0),
#'   method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42))
#' summary(fit)
#'
#' @seealso \code{\link{arma.cormat}}, \code{\link{poisson.marg}},
#'   \code{\link{zib.marg}}, \code{\link{zibb.marg}}, \code{\link{gctsc.opts}}
#' @export
#'
gctsc <- function(formula=NULL, data, marginal, cormat,
                  method = c("TMET", "GHK", "CE"),
                  c = 0.5, QMC = TRUE, pm = 30, start = NULL,
                 options = gctsc.opts()) {

  method <- match.arg(method)
  .validate_method(method, "gctsc")
  objs <- .validate_marg_cormat(marginal, cormat, "gctsc")
  marginal <- objs$marginal; cormat <- objs$cormat
  formula  <- .validate_formula_input(formula, marginal, "gctsc")
  .validate_options(method, QMC, options, "gctsc")
  
  des <- .build_design(formula, data, marginal, "gctsc")
  y <- des$y; x <- des$x
  validate_x_structure(x, marginal, "gctsc")
  check_x_nrow_matches_y(x, y, marginal, "gctsc")
  
  # hand off
  fit <- gctsc.fit(x = x, y = y, marginal = marginal, cormat = cormat,
                   method = method, c = c, QMC = QMC, pm = pm,
                   start = start, options = options)
  
  fit$call <- match.call(expand.dots = FALSE)
  fit$formula <- formula
  fit$terms <- des$terms
  fit$model <- des$model
  class(fit) <- "gctsc"
  fit
}


#' Fit a Gaussian Copula Time Series Model (Internal)
#'
#' Internal workhorse called by \code{\link{gctsc}}. Validates inputs, builds
#' starting values and bounds from the marginal and correlation structures, and
#' maximizes the approximate log–likelihood for the chosen method.
#'
#' @inheritParams gctsc
#' @param x Design matrix (non–ZI) or list of design matrices \code{list(mu = X_mu, pi0 = X_pi0)} (ZI).
#' @param y Numeric response vector of non–negative integer counts.
#' @return A list with estimates, log–likelihood, (optionally) Hessian, and diagnostics.
#' @keywords internal
#' @seealso \code{\link{gctsc}}
#' @noRd
gctsc.fit <- function(x = NULL, y, marginal, cormat,
                      method = "GHK", c = 0.5, QMC = TRUE,
                      start = NULL, pm = 30, options = gctsc.opts()) {

  objs <- .validate_marg_cormat(marginal, cormat, "gctsc.fit")
  marginal <- objs$marginal; cormat <- objs$cormat
  .validate_method(method, "gctsc.fit")
  .validate_options(method, QMC, options, "gctsc.fit")
  
  if (is.null(x)) {
    if (has_ZI(marginal)) {
      x <- list(mu = matrix(1, nrow = length(y), ncol = 1L),
                pi0 = matrix(1, nrow = length(y), ncol = 1L))
    } else {
      x <- matrix(1, nrow = length(y), ncol = 1L)
    }
  }
  validate_x_structure(x, marginal, "gctsc.fit")
  check_x_nrow_matches_y(x, y, marginal, "gctsc.fit")
  
  # Missing handling
  x_mat <- if (has_ZI(marginal)) do.call(cbind, x) else x
  not.na <- rowSums(is.na(cbind(y, x_mat))) == 0
  if (!any(not.na)) stop("gctsc.fit(): Have NA after combining y and x.", call. = FALSE)
  if (sum(!not.na) > 0) warning(sprintf("gctsc.fit(): dropping %d row(s) with NA.", sum(!not.na)))
  
  y <- as.matrix(y)[not.na, , drop = FALSE]
  x <- if (has_ZI(marginal)) {
    lapply(x, function(col) as.matrix(col)[not.na, , drop = FALSE])
  } else {
    as.matrix(x)[not.na, , drop = FALSE]
  }
  
  nbeta <- marginal$npar(x)
  ntau  <- cormat$npar
  
  # Starting values & bounds (always read template attrs)
  beta_tmpl <- marginal$start(y, x)
  tau_tmpl  <- cormat$start(y)
  lb <- c(attr(beta_tmpl, "lower") %||% rep(-Inf, length(beta_tmpl)),
          attr(tau_tmpl,  "lower") %||% rep(-Inf, length(tau_tmpl)))
  ub <- c(attr(beta_tmpl, "upper") %||% rep( Inf, length(beta_tmpl)),
          attr(tau_tmpl,  "upper") %||% rep( Inf, length(tau_tmpl)))
  
  if (is.null(start)) {
    init_eta <- c(beta_tmpl, tau_tmpl)
  } else {
    if (!is.numeric(start) || length(start) != (nbeta + ntau) || any(!is.finite(start)))
      stop("gctsc.fit(): 'start' must be numeric, finite, length nbeta + ntau.", call. = FALSE)
    init_eta <- start
  }
  
  # Optional AR/MA admissibility check at initial tau
  if (ntau > 0 && !is.null(cormat$p) && !is.null(cormat$q)) {
    p <- cormat$p; q <- cormat$q
    tau_init <- tail(init_eta, ntau)
    .check_arima_admissibility(tau_init, p, q, "gctsc.fit")
  }
  
  f <- structure(list(
    y = y, x = x, c = c, n = sum(not.na), method = method,
    marginal = marginal, cormat = cormat,
    ibeta = 1:nbeta, itau = (nbeta + 1):(nbeta + ntau),
    nbeta = nbeta, ntau = ntau, QMC = QMC, pm = pm,
    call = match.call(), init_eta = init_eta, coef = init_eta,
    lower = lb, upper = ub, options = options
  ), class = "gctsc")
  
  gctsc.estimate(f)
}



#' @keywords internal
#' @noRd
gctsc.estimate <- function(cf) {

  
  if (!is.null(cf$options$seed)) {
    # Temporarily set user-provided seed
    set.seed(cf$options$seed)
  }
  

  start <- cf$init_eta
  low <- cf$lower
  up <- cf$upper
  penalty <- -sqrt(.Machine$double.xmax)
  M <- cf$options$M
  log.lik <- build_loglik(cf,M, penalty)

  # saving/restoring the random seed (only for methods that need it)
  ans <- suppressWarnings(cf$options$opt(start, log.lik, low, up))
  

  eta <- ans$estimate
  names(eta) <- names(cf$coef)
  cf$coef <- eta
  cf$maximum <- ans$maximum
  cf$convergence <- ans$convergence

  # Store Hessian if available
  if (!is.null(ans$hessian) && is.matrix(ans$hessian) && all(is.finite(ans$hessian))) {
    cf$hessian <- ans$hessian
    if (all(eigen(cf$hessian, symmetric = TRUE, only.values = TRUE)$values > 0)) {
      vcov <- try(solve(cf$hessian), silent = TRUE)
      if (!inherits(vcov, "try-error")) {
        cf$se <- sqrt(diag(vcov))
      }
    } else {
      warning("Hessian is not positive definite. Standard errors not computed.")
      cf$se <- rep(NA, length(cf$coef))
    }
  } else {
    warning("Hessian not available from optimization.")
    cf$se <- rep(NA, length(cf$coef))
  }

  return(cf)
}








