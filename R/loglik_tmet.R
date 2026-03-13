#' Approximate Log-Likelihood via TMET
#'
#' Computes the approximate log-likelihood for a count time series model
#' based on a Gaussian and Student --t copula using the
#' \emph{Time Series Minimax Exponential Tilting (TMET)} method.
#'
#' TMET exploits the autoregressive moving-average (ARMA) structure
#' of the latent Gaussian or Scale mixture normal representation of
#' Student--t process to evaluate high-dimensional
#' multivariate normal/t rectangle probabilities via adaptive
#' importance sampling with an optimal tilting parameter.
#'
#' The implementation combines the Innovations Algorithm for exact
#' conditional mean and variance computation with exponential tilting,
#' resulting in a scalable and variance-efficient likelihood approximation.
#'
#' In this package, the latent dependence structure is parameterized
#' through an ARMA(\eqn{p,q}) process.
#'
#' @param lower Numeric vector of length \code{n} giving the lower
#'   bounds of the transformed latent variables.
#'
#' @param upper Numeric vector of length \code{n} giving the upper
#'   bounds of the transformed latent variables.
#'
#' @param tau Numeric vector of ARMA dependence parameters ordered as
#'   \code{c(phi_1, ..., phi_p, theta_1, ..., theta_q)}, where
#'   \eqn{\phi_i} are autoregressive (AR) coefficients and
#'   \eqn{\theta_j} are moving-average (MA) coefficients.
#'
#' @param od Integer vector \code{c(p, q)} specifying the AR and MA orders.
#'
#' @param M Positive integer specifying the number of Monte Carlo or
#'   quasi-Monte Carlo samples used in the simulation.
#'
#' @param QMC Logical; if \code{TRUE} (default), quasi-Monte Carlo
#'   integration is used. Otherwise, standard Monte Carlo sampling
#'   is applied.
#'
#' @param ret_llk Logical; if \code{TRUE} (default), returns the approximate
#' log-likelihood. If \code{FALSE}, internal diagnostic quantities
#' from the GHK simulator are returned. This option is primarily
#' intended for internal use and methodological research.
#'
#' @param df Degrees of freedom for the t copula. Must be greater than 2.
#'   Required only for \code{pmvt_ghk()}.
#' @param pm Integer specifying the number of past lags used when
#'   approximating an ARMA(\eqn{p,q}) process by an AR representation
#'   (required if \code{q > 0}).
#'
#' @return A numeric scalar giving the approximate log-likelihood.
#'   If \code{ret_llk = FALSE}, diagnostic output from the TMET
#'   sampler is returned (primarily for research use).
#' @seealso \code{\link{pmvn_ghk}}, \code{\link{pmvt_ghk}}
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
#' ## Gaussian copula example
#' mu <- 10
#' tau <- 0.2
#' arma_order <- c(1, 0)
#'
#' sim_data <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                         nsim = 1000, seed = 1)
#'
#' y <- sim_data$y
#' a <- qnorm(ppois(y - 1, lambda = mu))
#' b <- qnorm(ppois(y, lambda = mu))
#'
#' # Approximate log-likelihood using TMET
#' llk_tmet <- pmvn_tmet(lower = a, upper = b,
#'                       tau = tau, od = arma_order)
#' llk_tmet
#' 
#' ## Student--t copula example
#' df <- 8
#'
#' sim_data_t <- sim_poisson(mu = mu, tau = tau, arma_order = arma_order,
#'                           nsim = 500, family = "t", df = df, seed = 1)
#'
#' y_t <- sim_data_t$y
#' a_t <- qt(ppois(y_t - 1, lambda = mu), df = df)
#' b_t <- qt(ppois(y_t, lambda = mu), df = df)
#'
#' llk_t <- pmvt_tmet(lower = a_t, upper = b_t, tau = tau, od = arma_order,
#'                   M = 1000, df = df)
#'                   
#' @rdname pmv_tmet
#' @export
pmvn_tmet <- function(lower, upper, tau, od,
                      pm = 30, M = 1000,
                      QMC = TRUE, ret_llk = TRUE) {
  
  tmet_core(lower, upper, tau, od,
            family = "gaussian",
            pm = pm, M = M,
            QMC = QMC,
            ret_llk = ret_llk)
}

#' @rdname pmv_tmet 
#' @export
pmvt_tmet <- function(lower, upper, tau, od,
                      pm = 30, M = 1000,
                      QMC = TRUE, ret_llk = TRUE,
                      df) {
  
  tmet_core(lower, upper, tau, od,
            family = "t",
            pm = pm, M = M,
            QMC = QMC,
            ret_llk = ret_llk,
            df = df)
}

#' @keywords internal
#' @noRd
loglik_tmet <- function(ab, tau, od,
                        family,
                        pm = 30,
                        M = 1000,
                        QMC = TRUE,
                        ret_llk = TRUE,
                        df = NULL) {
  
  if (any(is.na(ab)))
    return(-1e20)
  
  tryCatch(
    tmet_core(ab[,1], ab[,2], tau, od,
              family = family,
              pm = pm,
              M = M,
              QMC = QMC,
              ret_llk = ret_llk,
              df = df),
    error = function(e) -1e20
  )
}


tmet_core <- function(lower, upper, tau, od,
                      family,
                      pm = 30,
                      M = 1000,
                      QMC = TRUE,
                      ret_llk = TRUE,
                      df = NULL) {
  
  if (length(tau) != sum(od))
    stop("Length of 'tau' must equal p+q.")
  
  if (all(od == 0))
    stop("ARMA(0,0) not supported.")
  
  if (any(upper < lower))
    stop("Invalid bounds: upper < lower.")
  
  n <- length(lower)
  p <- od[1]
  q <- od[2]
  
  pm <- if (q == 0) p else pm
  NN <- build_NN(n, pm)
  
  tmet_obj <- cond_mv_tmet(NN, tau, od)
  
  if (family == "gaussian") {
    
    # ----- Gaussian tilting solve -----
    z0 <- truncnorm::etruncnorm(lower, upper)
    z0_delta0 <- c(z0, rep(0, n))

    
    solv_delta <- stats::optim(
      z0_delta0,
      fn = function(x, ...) {
        ret <- grad_jacprod(x, ..., retProd = FALSE)
        0.5 * sum(ret$grad^2)
      },
      gr = function(x, ...) {
        ret <- grad_jacprod(x, ..., retProd = TRUE)
        ret$jac_grad
      },
      method = "L-BFGS-B",
      Condmv_Obj = tmet_obj,
      a = lower, b = upper,
      lower = c(lower, rep(-Inf, n)),
      upper = c(upper, rep(Inf, n)),
      control = list(maxit = 100000)
    )
    
    delta <- solv_delta$par[(n + 1):(2 * n)]
    
    exp_psi <- sample_mvn(
      tmet_obj, lower, upper,
      delta = delta,
      M = M, QMC = QMC,
      ret_llk = ret_llk,
      method = "TMET"
    )
    
    if (ret_llk) {
      exponent <- exp_psi[[2]]
      return(exponent + log(exp_psi[[1]]))
    } else {
      return(exp_psi$summary_stats)
    }
    
  } else {
    
    # ----- Student-t tilting solve -----
    trunc_expect <- truncnorm::etruncnorm(lower, upper)
    z0 <- c(trunc_expect, rep(0, n))
    z0[2 * n] <- 0
    z0[n] <- 1
    
    
    sol <- lm_sparse_solver(
      z0,
      Condmv_Obj = tmet_obj,
      a = lower, b = upper,
      nu = df,
      maxit = 500,
      tol = 1e-1,
      verbose = FALSE
    )
    
    z <- sol$x
    delta <- z[(n + 1):(2 * n)]
    eta <- delta[n]
    
    const <- log(2*pi)/2 - lgamma(df/2) -
      (df/2 - 1)*log(2) +
      TruncatedNormal::lnNpr(-eta, Inf) +
      0.5 * eta^2
    
    exp_psiT <- sample_mvt(
      tmet_obj, lower, upper,
      delta = delta,
      M = M, QMC = QMC,
      ret_llk = ret_llk,
      df = df,
      method = "TMET"
    )
    
    if (ret_llk) {
      exponent <- exp_psiT[[2]]
      log_est_prob <- exponent +log(mean(exp_psiT[[1]] * exp(exp_psiT[[2]] - exponent))) + const
      
      return(log_est_prob)
    } else {
      return(exp_psiT$summary_stats)
    }
  }
}


build_NN <- function(n, pm) {
  NN <- matrix(NA_integer_, nrow = n, ncol = pm + 1)
  NN[1, 1] <- 1
  
  if (n > 1) {
    row_idx <- 2:n
    col_idx <- 0:pm
    
    idx_mat <- outer(row_idx - 1, col_idx, FUN = function(i, j) i - j)
    valid_mask <- outer(row_idx - 1, col_idx, FUN = function(i, j) j <= i)
    idx_mat[!valid_mask] <- NA_integer_
    
    NN[2:n, ] <- idx_mat + 1
  }
  
  NN
}