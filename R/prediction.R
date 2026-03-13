#' @name predict.gctsc
#' @title One-Step-Ahead Predictive Distribution for Copula Count Time Series Models
#'
#' @description
#' Computes the one-step-ahead predictive distribution for a fitted
#' Gaussian or Student--t copula count time series model.
#'
#' The predictive probability mass function is evaluated using the
#' estimation method stored in the fitted object (TMET or GHK).
#' Summary statistics of the predictive distribution are returned,
#' and optional scoring rules are computed if the observed value is supplied.
#'
#' @param object A fitted model object of class \code{"gctsc"},
#'   as returned by \code{\link{gctsc}}.
#'
#' @param y_obs Optional non-negative integer giving the observed value
#'   at the prediction time. If supplied, the Continuous Ranked Probability
#'   Score (CRPS) and Logarithmic Score (LOGS) are computed.
#'
#' @param X_test Covariate values at the prediction time point.
#'   Can be provided as:
#'   \itemize{
#'     \item a named numeric vector whose names match the covariates
#'           used during model fitting, or
#'     \item a one-row \code{data.frame} or matrix with matching column names.
#'   }
#'   If the fitted model is intercept-only, \code{X_test} may be omitted,
#'   in which case the intercept is set to 1 automatically.
#'   An error is raised if more than one row is supplied.
#' @param ... Ignored. Included for S3 method compatibility.
#' 
#' @return A list containing:
#' \itemize{
#'   \item \code{mean}: Predictive mean.
#'   \item \code{median}: Predictive median.
#'   \item \code{mode}: Predictive mode.
#'   \item \code{variance}: Predictive variance.
#'   \item \code{p_y}: Predictive probability mass function over
#'         \code{0:y_max}.
#'   \item \code{lower}, \code{upper}: Bounds of the 95\% predictive interval.
#'   \item \code{CRPS}: Continuous Ranked Probability Score
#'         (if \code{y_obs} is provided).
#'   \item \code{LOGS}: Logarithmic Score
#'         (if \code{y_obs} is provided).
#' }
#'
#' @details
#' The predictive distribution integrates over the latent copula
#' dependence structure specified in the fitted model. For Gaussian
#' copulas, multivariate normal rectangle probabilities are evaluated;
#' for Student--t copulas, multivariate t rectangle probabilities are used.
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
#' @examples
#' # Simulate Poisson AR(1) data
#' set.seed(1)
#' y_sim <- sim_poisson(mu = 10, tau = 0.2,
#'                      arma_order = c(1, 0),
#'                      nsim = 1000,
#'                      family = "gaussian")$y
#'
#' # Fit Gaussian copula model
#' fit <- gctsc(
#'   y ~ 1,
#'   data = data.frame(y = y_sim),
#'   marginal = poisson.marg(lambda.lower = 0),
#'   cormat = arma.cormat(p = 1, q = 0),
#'   method = "CE",
#'   family = "gaussian",
#'   options = gctsc.opts(M = 1000, seed = 42)
#' )
#'
#' # One-step-ahead prediction
#' predict(fit)
#'
#' @seealso \code{\link{gctsc}}, \code{\link{arma.cormat}}
#' @export
#' @method predict gctsc
#'
predict.gctsc <- function(object, y_obs = NULL, X_test = NULL,...) {
  if (!inherits(object, "gctsc"))
    stop("predict.gctsc() is only for objects of class 'gctsc'.")
  
  bounds <- object$marginal$bounds
  y <- object$y
  x <- object$x
  seed <- object$options$seed
  family <- object$family
  method <- object$method
  # Determine y_max for predictive distribution
  k <- 3
  y_sd <- sd(y)
  y_max <- floor(max(y) + k * y_sd)
  
  # --- Standardize X_test ---
  if (!is.null(X_test)) {
    # Convert 1-row data.frame or matrix to named numeric vector
    if (is.data.frame(X_test) || (is.matrix(X_test) && nrow(X_test) == 1)) {
      if (nrow(X_test) != 1) {
        stop("'X_test' must have exactly one row for prediction.")
      }
      nm <- colnames(X_test)  # store names first
      X_test <- as.numeric(X_test)
      names(X_test) <- nm
    }
    if (!is.numeric(X_test) || is.null(names(X_test))) {
      stop("'X_test' must be a named numeric vector, 1-row matrix, or 1-row data.frame.")
    }
  }
  
  # --- Fill in defaults or add intercepts ---
  if (is.null(X_test)) {
    if (has_ZI(object$marginal)) {
      if (has_only_intercept(x$mu) && has_only_intercept(x$pi0)) {
        X_test <- c(
          setNames(rep(1, ncol(x$mu)),  colnames(x$mu)),
          setNames(rep(1, ncol(x$pi0)), colnames(x$pi0))
        )
      } else {
        stop("'X_test' must be provided as a named numeric vector matching covariate names.")
      }
    } else {
      if (has_only_intercept(x)) {
        X_test <- setNames(rep(1, ncol(x)), colnames(x))
      } else {
        stop("'X_test' must be provided as a named numeric vector matching covariate names.")
      }
    }
  } else {
    if (has_ZI(object$marginal)) {
      req_mu  <- colnames(x$mu)
      req_pi0 <- colnames(x$pi0)
      X_test <- add_intercept_if_needed(X_test, req_mu)
      X_test <- add_intercept_if_needed(X_test, req_pi0)
      
      if (!all(req_mu %in% names(X_test))) {
        stop("Missing mu covariates: ", paste(setdiff(req_mu, names(X_test)), collapse = ", "))
      }
      if (!all(req_pi0 %in% names(X_test))) {
        stop("Missing pi0 covariates: ", paste(setdiff(req_pi0, names(X_test)), collapse = ", "))
      }
    } else {
      req_mu <- colnames(x)
      X_test <- add_intercept_if_needed(X_test, req_mu)
      
      if (!all(req_mu %in% names(X_test))) {
        stop("Missing covariates: ", paste(setdiff(req_mu, names(X_test)), collapse = ", "))
      }
    }
  }
  
  # --- Construct X_rep ---
  if (has_ZI(object$marginal)) {
    X_mu  <- matrix(X_test[colnames(x$mu)],  nrow = 1)
    X_pi0 <- matrix(X_test[colnames(x$pi0)], nrow = 1)
    colnames(X_mu)  <- colnames(x$mu)
    colnames(X_pi0) <- colnames(x$pi0)
    
    X_rep <- list(
      mu  = X_mu[rep(1, y_max + 1), , drop = FALSE],
      pi0 = X_pi0[rep(1, y_max + 1), , drop = FALSE]
    )
  } else {
    X_test_mat <- matrix(X_test[colnames(x)], nrow = 1)
    colnames(X_test_mat) <- colnames(x)
    X_rep <- X_test_mat[rep(1, y_max + 1), , drop = FALSE]
  }
  

  # Bounds for observed y and prediction values
  ab   <- bounds(y, x, object$coef[object$ibeta], family = family, df= object$df)
  ab_p <- bounds(0:y_max, X_rep, object$coef[object$ibeta], family = family, df= object$df)
  
  # Predictive distribution setup
  pred_input <- list(
    a     = ab[, 1],
    b     = ab[, 2],
    ap    = ab_p[, 1],
    bp    = ab_p[, 2],
    y_max = y_max,
    M     = object$options$M,
    QMC   = object$QMC
  )
  
  od  <- object$cormat$od
  tau <- object$coef[object$itau]
  set.seed(seed)
  if (family =="gaussian"){
    if (method == "TMET") {
      p_y <- pred_tmet_mvn(pred_input, tau = tau, od = od)$p_y
    } else {
      p_y <- pred_ghk_mvn(pred_input, tau = tau, od = od)$p_y
    }
  } else{
    pred_input$df <- object$df
    p_y <- pred_mvt(args = pred_input, tau = tau, od = od)$p_y
    }
  
  
  # Summary stats
  y_vals      <- 0:y_max
  cdf         <- cumsum(p_y)
  mean_pred   <- sum(y_vals * p_y)
  var_pred    <- sum((y_vals^2) * p_y) - mean_pred^2
  median_pred <- y_vals[which(cdf >= 0.5)[1]]
  mode_pred   <- y_vals[which.max(p_y)]
  
  alpha <- 0.05
  cdf_shifted <- c(0, cdf[-length(cdf)])
  
  lower <- ifelse(
    any(cdf_shifted <= alpha / 2),
    y_vals[max(which(cdf_shifted <= alpha / 2))],
    y_vals[1]
  )
  upper <- ifelse(
    any(cdf >= 1 - alpha / 2),
    y_vals[min(which(cdf >= 1 - alpha / 2))],
    y_vals[length(y_vals)]
  )
  
  out <- list(
    mean    = mean_pred,
    median  = median_pred,
    mode    = mode_pred,
    variance = var_pred,
    p_y     = p_y,
    lower   = lower,
    upper   = upper
  )
  
  if (!is.null(y_obs)) {
    indicator <- as.numeric(y_vals >= y_obs)
    out$CRPS  <- sum((cdf - indicator)^2)
    out$LOGS  <- -log(p_y[y_obs + 1])
  }
  
  return(out)
}


#' @keywords internal
#' @noRd
pred_tmet_mvn <- function(args, tau, od){
  if (length(tau) != sum(od))
    stop("Length of 'tau' must match ARMA order")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(is.na(args$a)) || any(is.nan(args$a))) return(list(Ey = NA, VEy = NA))
  
  
  .p <- od[1]
  .q <- od[2]
  iar <- if (.p) 1:.p else NULL
  ima <- if (.q) (.p + 1):(.p + .q) else NULL
  phi <- tau[iar]
  theta <- tau[ima]
  
  if(.p==0){
    p<- 1
    phi <-0
  } else {p<-.p}
  if(.q==0){
    q<-1
    theta<-0
  } else{q <- .q}
  
  m = max(p,q)
  pm <- if (.q == 0) p else 30
  Tau <- list(phi = phi, theta = theta)
  n <- length(args$a)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))
  
  
  
  
  ## NNarray: causal lag indices  --------------------------------
  
  NN <- matrix(NA_integer_, nrow = n, ncol = pm + 1)
  NN[1, 1] <- 1
  
  row_idx <- 2:n
  col_idx <- 0:pm
  
  # Create a matrix of row indices and subtract column offsets
  idx_mat <- outer(row_idx - 1, col_idx, FUN = function(i, j) i - j)
  
  # Mask values where j > i (i.e., invalid entries)
  valid_mask <- outer(row_idx - 1, col_idx, FUN = function(i, j) j <= i)
  idx_mat[!valid_mask] <- NA_integer_
  
  NN[2:n, ] <- idx_mat +1
  
  
  # compute conditional variance and BLUP coefficient mt B
  tmet_obj <- cond_mv_tmet(NN,tau,od)
  cond_sd <- sqrt(tmet_obj$cond_var)
  Theta <-  rbind(tmet_obj$Theta)
 
  
  # find tilting parameter delta -----------------------------------
  lower = args$a
  upper = args$b
  z0 <- truncnorm::etruncnorm(lower, upper)
  z0_delta0 <- c(z0, rep(0, n))
  
  solv_delta <- stats::optim(
    z0_delta0,
    fn = function(x, ...) {
      ret <- grad_jacprod(x, ...,retProd = FALSE)
      0.5*sum((ret$grad)^2)
    },
    gr = function(x, ...) {
      ret <- grad_jacprod(x, ..., retProd = TRUE)
      ret$jac_grad
    },
    method = "L-BFGS-B",
    Condmv_Obj = tmet_obj,
    a = lower, b = upper, 
    lower = c(lower, rep(-Inf, n)), upper = c(upper, rep(Inf, n)),
    control = list(maxit = 500)
  )
  
  if (any(solv_delta$par[1:n] < lower) ||
      any(solv_delta$par[1:n] > upper)) {
    warning("Optimal x is outside the integration region during minmax tilting\n")
  }
  
  delta<- solv_delta$par[(n + 1):(2 * n)]
  
  
  model <- list(
    phi = phi, theta_r = theta_r,
    n = n, p = p, q = q, m = m,
    sigma2 = sigma2, gamma = gamma, delta =delta, Theta = Theta, condSd= cond_sd,
    v = tmet_obj$cond_var
  )
  
  result <- predmvn_tmet(args, model)
  
}


#' @keywords internal
#' @noRd
pred_ghk_mvn <- function(args, tau, od) {
  if (length(tau) != sum(od))
    stop("Length of 'tau' must match ARMA order")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(is.na(args$a)) || any(is.nan(args$a))) return(list(Ey = NA, VEy = NA))
  
  .p <- od[1]
  .q <- od[2]
  iar <- if (.p) 1:.p else NULL
  ima <- if (.q) (.p + 1):(.p + .q) else NULL
  phi <- tau[iar]
  theta <- tau[ima]
  
  if(.p==0){
    p<- 1
    phi <-0
  } else {p<-.p}
  if(.q==0){
    q<-1
    theta<-0
  } else{q <- .q}
  
  m = max(p,q)
  
  Tau <- list(phi = phi, theta = theta)
  n <- length(args$a)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))
  
  model <- list(
    phi = phi, theta_r = theta_r,
    n = n, p = p, q = q, m = m,
    sigma2 = sigma2, gamma = gamma
  )
  
  result <- predmvn_ghk(args, model)
}



#' @keywords internal
#' @noRd
pred_mvt <- function(args, tau, od) {
  if (length(tau) != sum(od))
    stop("Length of 'tau' must match ARMA order")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(is.na(args$a)) || any(is.nan(args$a))) return(list(Ey = NA, VEy = NA))
  
  .p <- od[1]
  .q <- od[2]
  iar <- if (.p) 1:.p else NULL
  ima <- if (.q) (.p + 1):(.p + .q) else NULL
  phi <- tau[iar]
  theta <- tau[ima]
  
  if(.p==0){
    p<- 1
    phi <-0
  } else {p<-.p}
  if(.q==0){
    q<-1
    theta<-0
  } else{q <- .q}
  
  m = max(p,q)
  
  Tau <- list(phi = phi, theta = theta)
  n <- length(args$a)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)
  gamma <- aacvf(Tau, n - 1)
  theta_r <- c(1, theta, numeric(n))
  
  model <- list(
    phi = phi, theta_r = theta_r,
    n = n, p = p, q = q, m = m,
    sigma2 = sigma2, gamma = gamma
  )
  
  result <- predmvt(args, model)
}

