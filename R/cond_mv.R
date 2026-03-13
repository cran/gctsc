#' @keywords internal
#' @noRd
cond_mv_ghk <- function(n, tau, od) {
  .p <- od[1]
  .q <- od[2]
  iar <- if ( .p ) 1:.p else NULL
  ima <- if ( .q ) (.p+1):(.p+.q) else NULL
  phi <- tau[iar]
  theta<-tau[ima]
  
  # Use default (p = 1, q = 1) if either AR or MA order is zero
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
  
  list(
    cond_var = cond_var,
    Theta = Theta,
    m =m,
    phi = phi,
    p = .p,
    q = .q
  )
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



