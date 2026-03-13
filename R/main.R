#' Fit a Copula-Based Count Time Series Model
#'
#' Fits a Gaussian or Student--t copula model to univariate count
#' time series with discrete marginals (Poisson, negative binomial,
#' binomial/beta–binomial, and zero–inflated variants) and latent
#' dependence specified via ARMA correlation structures.
#'
#' The high-dimensional rectangle probability defining the copula
#' likelihood is approximated using one of:
#' \itemize{
#'   \item \strong{TMET} (Time Series Minimax Exponential Tilting),
#'   \item \strong{GHK} (Geweke--Hajivassiliou--Keane simulation), or
#'   \item \strong{CE} (Continuous Extension).
#' }
#'
#' The interface mirrors \code{glm()}. Zero–inflated marginals accept a
#' named list of formulas, e.g.,
#' \code{list(mu = y ~ x, pi0 = ~ z)}. Non–zero–inflated marginals accept
#' a single formula or \code{list(mu = ...)}.
#'
#' @param formula A formula (e.g., \code{y ~ x1 + x2}) or, for
#'   zero–inflated marginals, a named list
#'   \code{list(mu = ..., pi0 = ...)}. 
#'
#' @param data A data frame containing the response and covariates
#'   referenced in the formula(s).
#'
#' @param marginal A marginal model object such as
#'   \code{\link{poisson.marg}},
#'   \code{\link{negbin.marg}},
#'   \code{\link{zib.marg}}, or
#'   \code{\link{zibb.marg}}; must inherit class
#'   \code{"marginal.gctsc"}.
#'
#' @param cormat A correlation structure such as
#'   \code{\link{arma.cormat}}; must inherit class
#'   \code{"cormat.gctsc"}.
#'
#' @param method One of \code{"TMET"}, \code{"GHK"}, or \code{"CE"}.
#'
#' @param c Smoothing constant used by CE only (ignored otherwise).
#'   Default is \code{0.5}.
#'
#' @param QMC Logical; if \code{TRUE}, quasi–Monte Carlo integration
#'   is used for simulation–based methods.
#'
#' @param pm Integer specifying the truncated AR order used when
#'   approximating ARMA(\eqn{p,q}) by an AR representation
#'   (TMET only; relevant when \eqn{q > 0}). Default is 30
#'
#' @param start Optional numeric vector of starting values (marginal
#'   parameters followed by dependence parameters). If \code{NULL},
#'   sensible starting values and bounds are constructed from
#'   \code{marginal} and \code{cormat}.
#'
#' @param options Optional list of tuning and optimization controls.
#'   If \code{NULL}, defaults from \code{gctsc.opts()} are used
#'   (e.g., \code{M = 1000}, randomized seed). Any supplied fields
#'   override the defaults.
#' 
#' @param family Copula family. One of \code{"gaussian"}
#'   or \code{"t"}.
#'
#' @param df Degrees of freedom for the Student--t copula.
#'   Must be greater than 2. Required when \code{family = "t"}.
#'   Ignored for the Gaussian copula.
#'   
#' @details
#' \strong{Formulas.}
#' For zero–inflated marginals, if neither \code{mu} nor \code{pi0}
#' is supplied, both default to intercept–only models
#' (\code{mu ~ 1}, \code{pi0 ~ 1}). If \code{mu} is supplied but
#' \code{pi0} is missing, \code{pi0 ~ 1} is used.
#'
#' \strong{Dependence.}
#' The ARMA parameters are encoded in \code{cormat}. Models must
#' satisfy stationarity and invertibility conditions.
#' ARMA(0,0) is not supported.
#'
#' \strong{Method-specific notes.}
#' CE ignores \code{QMC} and \code{options$M}.
#' GHK and TMET require \code{options$M} to be a positive integer.
#' TMET additionally uses \code{pm} when \eqn{q > 0}.
#'
#' @return An object of class \code{"gctsc"} containing, among others:
#' \itemize{
#'   \item \code{coef}: parameter estimates,
#'   \item \code{maximum}: approximate log–likelihood at the optimum,
#'   \item \code{se}: standard errors when available,
#'   \item \code{terms}, \code{model}, \code{call}: model metadata.
#' }
#'
#' @references
#' Nguyen, Q. N., & De Oliveira, V. (2026). Likelihood Inference in Gaussian Copula Models for Count Time Series
#' via Minimax Exponential Tilting
#' \emph{Journal of Computational Statistics and Data Analysis}.
#' 
#' Nguyen, Q. N., & De Oliveira, V. (2026).
#' Scalable Likelihood Inference for Student--\eqn{t} Copula Count Time Series.
#' Manuscript in preparation.
#' 
#' Nguyen, Q. N., & De Oliveira, V. (2026).
#' Approximating Gaussian copula models for count time series:
#' Connecting the distributional transform and a continuous extension.
#' \emph{Journal of Applied Statistics}.
#' 
#' @examples
#' ## Example 1: Gaussian copula, Poisson marginal, AR(1)
#' set.seed(42)
#' n <- 500
#' sim_dat <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1, 0),
#'                        nsim = n, family = "gaussian")
#'
#' dat <- data.frame(y = sim_dat$y)
#'
#' fit_gauss <- gctsc(
#'   y ~ 1,
#'   data = dat,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0), family = "gaussian",
#'   method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#' summary(fit_gauss)
#'
#' ## Example 2: Student--t copula
#' sim_dat_t <- sim_poisson(mu = 10, tau = 0.3, arma_order = c(1, 0),
#'                          nsim = 500, family = "t", df = 10)
#'
#' dat_t <- data.frame(y = sim_dat_t$y)
#'
#' fit_t <- gctsc(
#'   y ~ 1,
#'   data = dat_t,
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0), family ="t",
#'   df= 10, method = "CE",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#' summary(fit_t)
#'
#' @seealso \code{\link{arma.cormat}}, \code{\link{poisson.marg}},
#'   \code{\link{zib.marg}}, \code{\link{zibb.marg}}, \code{\link{gctsc.opts}}
#' @export
#'
gctsc <- function(formula=NULL, data, marginal, cormat,
                  method = c("TMET", "GHK", "CE"),
                  c = 0.5, QMC = TRUE, pm = 30, start = NULL,
                  family =c("t","gaussian"),df=10,
                  options = gctsc.opts()) {

  method <- match.arg(method)
  family <- match.arg(family)
  
  ## ---- t copula checks ----
  if (family == "t") {
    if (is.null(df))
      stop("For a Student-t copula, 'df' must be provided.")
    
    if (!is.numeric(df) || length(df) != 1)
      stop("'df' must be a single numeric value.")
    
    if (df <= 2)
      stop("'df' must be greater than 2 for the t copula.")
  }
  
  
  .validate_method(method, "gctsc")
  objs <- .validate_marg_cormat(marginal, cormat, "gctsc")
  marginal <- objs$marginal; cormat <- objs$cormat
  formula  <- .validate_formula_input(formula, marginal, "gctsc")
  .validate_options(method, QMC, options, "gctsc")
  
  des <- .build_design(formula, data, marginal, "gctsc")
  y <- des$y; x <- des$x
  validate_x_structure(x, marginal, "gctsc")
  check_x_nrow_matches_y(x, y, marginal, "gctsc")
  
  # passing family to marginal
  
  marginal$family <- family
  marginal$df <- df
  
  
  # hand off
  fit <- gctsc.fit(x = x, y = y, marginal = marginal, cormat = cormat,
                   method = method, c = c, QMC = QMC, pm = pm,df=df,
                   start = start, options = options, family = family)
  
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
                      method, c = 0.5, QMC = TRUE, df=NULL,
                      start = NULL, pm = 30, options = gctsc.opts(),
                      family = c("gaussian","t")) {

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
    lower = lb, upper = ub, options = options,family =family, df=df
  ), class = "gctsc")
  
  gctsc.estimate(f)
}





#' Set Options for Gaussian and Student t Copula Time Series Model
#'
#' Creates a control list for simulation and likelihood approximation in the
#' Gaussian ad Student t copula model, including the random seed and Monte Carlo settings.
#'
#' @param seed recommended to provde a seed. Integer specifying the random seed used for Monte Carlo
#' simulation during likelihood evaluation. Fixing the seed ensures
#' reproducible optimization results and stable maximum likelihood estimates.
#' @param M Integer. Number of Monte Carlo samples used in the likelihood approximation (default: 1000).
#' @param ... Ignored. Included for S3 method compatibility.
#'
#'
#' @return
#' A list with components:
#'   \item{\code{seed}}{ Integer. The random seed used.}
#'   \item{\code{M}}{ Integer. Number of Monte Carlo samples.}
#'   \item{\code{opt}}{ A function used internally by \code{gctsc()} to
#'         perform optimization of the approximate log-likelihood.}

#' 
#' @export
gctsc.opts <- function(seed=NULL, M=1000, ...) {
  
  
  control <- list(...)
  opt <- function(start, llk_fn, lower, upper) {
    fn.opt <- function(x)  {if( any(x <= lower | x >= upper) ) (1e10)
      # Compute the log-likelihood sum
      # else {-llk_fn(x)}
      else {
        val <- -llk_fn(x)
        if (!is.finite(val)) return(1e10)
        return(val)
      }
    }
    ans <- optim(start, fn.opt, method= "BFGS", hessian = TRUE)
    if(ans$convergence) warning(paste("optim exits with code",ans$convergence))
    list(
      estimate = ans$par,
      maximum = ans$value,
      convergence = ans$convergence,
      hessian = ans$hessian
    )
  }
  list(seed=seed,M=M,opt=opt)
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
  ans <- if (cf$method != "CE") {
    preserve_seed({
      suppressWarnings(cf$options$opt(start, log.lik, low, up))
    })
  } else {
    suppressWarnings(cf$options$opt(start, log.lik, low, up))
  }
  
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


