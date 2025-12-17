#' TMET Log-Likelihood Approximation
#'
#' Computes the approximate log-likelihood for a count time series model using the
#' \emph{Time Series Minimax Exponential Tilting (TMET)} method. This function uses
#' the autoregressive movin average structure of the underlying Gaussian copula to compute
#' multivariate normal rectangle probabilities via adaptive importance sampling
#' with an optimal tilting parameter.
#'
#' The implementation combines the Innovations Algorithm for exact conditional
#' mean and variance computation with exponential tilting to estimate
#' the multivariate normal probability over rectangular truncation sets.
#'
#' @param lower A numeric vector of lower truncation bounds.
#' @param upper A numeric vector of upper truncation bounds.
#' @param tau A numeric vector of ARMA dependence parameters.
#' @param od Integer vector \code{c(p, q)} specifying the AR and MA orders.
#' @param pm Integer. Number of past lags used for approximating ARMA(p,q) to AR  (required if \code{q > 0}).
#' @param M Integer. Number of Monte Carlo or QMC samples.
#' @param QMC Logical. If \code{TRUE}, use quasi-Monte Carlo integration.
#' @param ret_llk Logical. Default is \code{TRUE} to return log-likelihood; otherwise, return diagnostic output.
#'
#' @return If \code{ret_llk = TRUE}, a numeric scalar (log-likelihood); else a list containing diagnostic statistics.
#'
#' @examples
#' # Simulate Poisson AR(1) data
#' sim_data <- sim_poisson(mu =10, tau=0.2, arma_order=c(1,0), nsim = 1000, seed = 1)
#' mu=10
#' tau=0.2
#' arma_order=c(1,0)
#' sim_data <- sim_poisson(mu =mu, tau=tau, arma_order=arma_order, nsim = 1000, seed = 1)
#' y <- sim_data$y
#'
#' # Compute latent bounds for CE method
#' a <- qnorm(ppois(y - 1, lambda = mu))  # lower bound
#' b <- qnorm(ppois(y, lambda = mu))      # upper bound
#'
#' # Approximate log-likelihood with CE method
#' llk_tmet <- pmvn_tmet(lower = a, upper = b, tau = 0.2, od = c(1,0))
#' print(llk_tmet)
#'
#' @export
pmvn_tmet <- function(lower, upper, tau, od, pm = 30,M = 1000, QMC = TRUE, ret_llk=TRUE){
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported. Please specify at least one AR or MA term.")
  }

  n <-length(lower)
  p <- od[1]
  q <- od[2]

  pm <- if (q == 0) p else pm
  if (any(upper < lower)) {
    stop("Invalid MVN probability: some upper bounds < lower bounds at index ",
         which(upper < lower)[1])
  }


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


  # find tilting parameter delta -----------------------------------
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


  # compute MVN probs and est error ---------------------------------
  exp_psi <- sample_latent(tmet_obj, lower, upper, delta = delta, M = M,QMC=QMC, ret_llk = ret_llk,
                           method = "TMET")


  if(ret_llk){
    exponent <- exp_psi[[2]]
    log_pmvn <- exponent +log(exp_psi[[1]])
    return (log_pmvn)
    # lw <- exp_psi[[3]]
    # delta <- delta
    return(list("weights" = lw, "delta"=delta))
  }else{
    return(exp_psi$summary_stats)
  }
}

#' @keywords internal
#' @noRd
loglik_tmet<- function(ab, tau, M = 1000, od, pm = 30, QMC=TRUE, ret_llk = TRUE) {
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported. Please specify at least one AR or MA term.")
  }


  if (any(is.na(ab)) || any(is.nan(ab))) return(-1e20)
  a <- ab[,1]
  b <- ab[,2]
  p <- od[1]
  q <- od[2]

  pm <- if (q == 0) p else pm

  result <- tryCatch({
    pmvn_tmet(a, b, tau, od, pm,M,
                     QMC = QMC, ret_llk = ret_llk)
  }, error = function(e) {
    message("TMET failed: ", e$message)
    return(-1e20)
  })

  return(result)
}

#' @keywords internal
#' @noRd
cond_mv_tmet <- function(NN, tau, od) {
  n <- nrow(NN)
  if (!all(NN[, 1] == 1:n)) stop("Unexpected NN: first col is not 1:n\n")
  pm <- ncol(NN) - 1
  .p <- od[1]
  .q <- od[2]
  iar <- if ( .p ) 1:.p else NULL
  ima <- if ( .q ) (.p+1):(.p+.q) else NULL
  phi <- tau[iar]
  theta<-tau[ima]

  ### innovation algorithm assume p>=1, q>=1
  if(.p==0){
    p<- 1
    phi <-0
  } else {p<-.p}
  if(.q==0){
    q<-1
    theta<-0
  } else{q <- .q}

  m = max(p,q)
  Tau <- list(phi=phi, theta=theta)
  sigma2 <- 1 / sum(ma.inf(Tau)^2)

  gamma <- aacvf(Tau, n - 1)

  theta_r <- c(1, theta, numeric(n))

  model <- list(phi = phi, theta_r = theta_r, p = p, q = q, m = m,
                sigma2 = sigma2, gamma = gamma, n = n)
  out_cpp <- compute_cond_var(gamma, model)
  v <- out_cpp$v
  cond_var <- v*sigma2

  ##### theta_tk is saved to compute the conditional mean
  Theta <- out_cpp$Theta
  arma_blup_coef = - (ar.inf(Tau,pm)[-1])

  cond_mean_coeff <- matrix(0, n, pm)
  # first loop: i from 2:(pm-1)
  if(pm > 2){
    DL <- durbin_levinson(gamma, pm,cond_var[-1])
    for (i in 2:(pm - 1)) {
      cond_mean_coeff[i, 1:(i - 1)] <- DL[(i-1),1:(i - 1)]
    }
  }

  # Second loop: i from (pm - 1) to n
  for (i in max(2, pm):n) {
    cond_mean_coeff[i, 1:pm] <- arma_blup_coef[1:pm]
  }


  B <- sparse_B(NN, cond_mean_coeff, n, pm)

  list(
    cond_var = cond_var,
    cond_mean_coeff = cond_mean_coeff,
    B = B,
    NN = NN,
    Theta = Theta,
    phi = phi,
    p = .p,
    q = .q,
    m = m
  )

}



#' @keywords internal
#' @noRd
sparse_B <- function(NNarray, cond_mean_coeff, n, pm) {
  nnz_per_row <- pmin(pm + 1, 1:n)
  total_nnz <- sum(nnz_per_row)

  B_row_inds <- integer(total_nnz)
  B_col_inds <- integer(total_nnz)
  B_vals     <- numeric(total_nnz)

  pos <- 1
  for (i in 1:n) {
    k <- nnz_per_row[i]
    idx <- pos:(pos + k - 1)

    B_row_inds[idx] <- i
    B_col_inds[idx] <- NNarray[i, 1:k]

    if (k > 1) {
      B_vals[idx[-1]] <- cond_mean_coeff[i, seq_len(k - 1)]
    }

    pos <- pos + k
  }

  Matrix::sparseMatrix(i = B_row_inds, j = B_col_inds, x = B_vals, dims = c(n, n))
}







