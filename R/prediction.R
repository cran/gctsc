#' @name predict.gctsc
#' @aliases predict.gctsc predict
#' @title Predictive Distribution and Scoring for Gaussian Copula Time Series Models
#'
#' @description
#' Computes the one-step-ahead predictive distribution for a fitted Gaussian copula
#' time series model, including summary statistics (mean, median, mode, variance)
#' and optional scoring rules (CRPS and LOGS) if an observed value is supplied.
#'
#' @usage
#' \method{predict}{gctsc}(object, ..., method = "GHK",
#'                         y_obs = NULL, X_test = NULL)
#'
#' @param object A fitted model object of class \code{gctsc}.
#' @param ... Not used. Included for S3 method consistency.
#' @param method Character string specifying the prediction method:
#'   \code{"TMET"} or \code{"GHK"} (default: \code{"GHK"}).
#' @param y_obs Optional observed value used to compute CRPS and LOGS.
#' @param X_test Optional covariate information for prediction. Can be a named
#'   numeric vector, or a 1-row matrix/data.frame with column names matching the
#'   model covariates.
#'
#' @return 
#' \item{mean}{Predictive mean of \eqn{Y_{t+1}}.}
#' \item{median}{Predictive median of \eqn{Y_{t+1}}.}
#' \item{mode}{Predictive mode of \eqn{Y_{t+1}}.}
#' \item{variance}{Predictive variance of \eqn{Y_{t+1}}.}
#' \item{CRPS}{Continuous Ranked Probability Score (if \code{y_obs} is provided).}
#' \item{LOGS}{Logarithmic score (if \code{y_obs} is provided).}
#' \item{p_y}{Predictive pmf over \eqn{(0,\dots,y_{\max})}.}
#' \item{lower, upper}{ 95\eqn{%} predictive interval bounds.}
#' @seealso \code{\link{gctsc}}, \code{\link{arma.cormat}}
#' @export
#' @method predict gctsc
predict.gctsc <- function(object,
                          ...,
                          method = "GHK",
                          y_obs = NULL,
                          X_test = NULL) {
  
  if (!inherits(object, "gctsc"))
    stop("predict.gctsc() is only for objects of class 'gctsc'.")
  
  ## --------- Extract model components -----------
  bounds <- object$marginal$bounds
  y      <- object$y
  x      <- object$x
  seed   <- object$options$seed
  
  ## --------- Compute y_max -----------
  k     <- 3
  y_sd  <- sd(y)
  y_max <- floor(max(y) + k * y_sd)
  
  ## --------- Process X_test ------------
  # (your logic preserved exactly)
  if (!is.null(X_test)) {
    if (is.data.frame(X_test) || (is.matrix(X_test) && nrow(X_test) == 1)) {
      nm <- colnames(X_test)
      X_test <- as.numeric(X_test)
      names(X_test) <- nm
    }
    if (!is.numeric(X_test) || is.null(names(X_test))) {
      stop("'X_test' must be a named numeric vector, 1-row matrix, or 1-row data.frame.")
    }
  }
  
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
  
  ## --------- Construct X_rep -----------
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
  
  ## --------- Build predictive inputs ----------
  ab   <- bounds(y, x, object$coef[object$ibeta])
  ab_p <- bounds(0:y_max, X_rep, object$coef[object$ibeta])
  
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
  
  ## --------- Compute predictive PMF ----------
  p_y <- if (method == "TMET") {
    pred_tmet(pred_input, tau = tau, od = od)$p_y
  } else {
    pred_ghk(pred_input, tau = tau, od = od)$p_y
  }
  
  ## --------- Summary statistics ----------
  y_vals      <- 0:y_max
  cdf         <- cumsum(p_y)
  mean_pred   <- sum(y_vals * p_y)
  var_pred    <- sum(y_vals^2 * p_y) - mean_pred^2
  median_pred <- y_vals[which(cdf >= 0.5)[1]]
  mode_pred   <- y_vals[which.max(p_y)]
  
  alpha <- 0.05
  cdf_shifted <- c(0, cdf[-length(cdf)])
  
  lower <- ifelse(any(cdf_shifted <= alpha / 2),
                  y_vals[max(which(cdf_shifted <= alpha / 2))],
                  y_vals[1])
  upper <- ifelse(any(cdf >= 1 - alpha / 2),
                  y_vals[min(which(cdf >= 1 - alpha / 2))],
                  y_vals[length(y_vals)])
  
  out <- list(
    mean     = mean_pred,
    median   = median_pred,
    mode     = mode_pred,
    variance = var_pred,
    p_y      = p_y,
    lower    = lower,
    upper    = upper
  )
  
  ## --------- Add CRPS and LOGS if y_obs provided ----------
  if (!is.null(y_obs)) {
    indicator <- as.numeric(y_vals >= y_obs)
    out$CRPS  <- sum((cdf - indicator)^2)
    out$LOGS  <- -log(p_y[y_obs + 1])
  }
  
  return(out)
}


#' @keywords internal
#' @noRd
pred_tmet <- function(args, tau, od){
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
  cond_var <- tmet_obj$cond_var

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
    sigma2 = sigma2, gamma = gamma, delta =delta, Theta = Theta, condSd= cond_sd,v = tmet_obj$cond_var
  )

  result <- predmvn_tmet(args, model)

}

#' @keywords internal
#' @noRd
pred_ghk <- function(args, tau, od) {
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



