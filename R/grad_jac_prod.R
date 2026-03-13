# Compute gradient --------------------------------------------

#' @keywords internal
#' @noRd
grad_jacprod<- function(z_delta, Condmv_Obj, a, b, retProd = TRUE) {
  n <- length(a)
  z <- z_delta[1:n]
  delta <- z_delta[(n + 1):(2 * n)]
  D <- sqrt(Condmv_Obj$cond_var)
  B <- Condmv_Obj$B
  mu_c <- as.vector(B %*% z)
  a_t_shift <- (a - mu_c) / D - delta
  b_t_shift <- (b - mu_c) / D - delta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_t_shift, b_t_shift)
  pl <- exp(-0.5 * a_t_shift^2 - log_diff_cdf)  / sqrt(2 * pi)
  pu <- exp(-0.5 * b_t_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  ld <- pl - pu
  
  # compute grad ------------------------------------------------
  dpsi_dz <- as.vector(Matrix::t(B) %*%
                         (delta / D + ld / D)) - delta / D
  dpsi_ddelta <- delta - (z - mu_c) / D + ld
  
  # build return list ----------------------------------------------
  rslt <- list(
    grad = c(dpsi_dz, dpsi_ddelta)
  )
  if (retProd) {
    
    a_t_shift[is.infinite(a_t_shift)] <- 0
    b_t_shift[is.infinite(b_t_shift)] <- 0
    LD <- (-ld^2) + a_t_shift * pl - b_t_shift * pu
    
    H11_dpsi_dz <- as.vector(Matrix::t(LD / D / D * as.vector(B %*% dpsi_dz)) %*%B)
    H12_dpsi_ddelta <- as.vector(Matrix::t((dpsi_ddelta + LD * dpsi_ddelta) / D) %*% B) - dpsi_ddelta / D
    H21_dpsi_dz <-  as.vector(B %*% dpsi_dz) / D * (1 + LD) - dpsi_dz / D
    H22_dpsi_ddelta  <- ((1 + LD) * dpsi_ddelta)
    rslt$jac_grad <- c(
      H11_dpsi_dz + H12_dpsi_ddelta,
      H21_dpsi_dz + H22_dpsi_ddelta
    )
  }
  return(rslt)
}



grad_jac_psiT  <- function(z_delta, Condmv_Obj, a, b, nu, deriv = c("grad","jac","both")) {
  deriv <- match.arg(deriv)
  
  n <- length(a)
  z <- delta <- numeric(n)
  z[-n] <- z_delta[1:(n - 1)]
  w <- z_delta[n]
  delta[-n] <- z_delta[(n + 1):(2 * n - 1)]
  kap <- z_delta[2 * n]
  
  a <- a / sqrt(nu)
  b <- b / sqrt(nu)
  B <- Condmv_Obj$B
  D <- sqrt(Condmv_Obj$cond_var)
  
  mu_c <- as.vector(B %*% z)
  a_tilde_shift <- (a * w - mu_c) / D - delta
  b_tilde_shift <- (b * w - mu_c) / D - delta
  
  
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  ld <- pl - pu
  
  
  ## gradient ------------------------------------------------------
  if (deriv %in% c("grad","both")) {
    dpsi_dz <- as.vector(Matrix::t(B) %*% (delta / D + ld / D)) - delta / D
    dpsi_ddelta <- delta - (z - mu_c) / D + ld
    dpsi_dw <- (nu - 1) / w - kap + sum((b/D) * pu - (a/D) * pl)
    dpsi_dk <- kap - w + exp(stats::dnorm(kap, log = TRUE) -
                               TruncatedNormal::lnNpr(-kap, Inf))
    grad <- c(dpsi_dz[-n], dpsi_dw, dpsi_ddelta[-n], dpsi_dk)
  } else {
    grad <- NULL
  }
  
  ## Jacobian ------------------------------------------------------
  if (deriv %in% c("jac","both")) {
    dld_dw <- (-(a/D) * a_tilde_shift * pl) + (b/D) * b_tilde_shift * pu +
      (a/D) * pl^2 + (b/D) * pu^2 -
      ((a + b)/D) * (pl * pu)
    
    LD <- (-ld^2) + a_tilde_shift * pl - b_tilde_shift * pu
    
    # sparse helpers
    LD_diag   <- Matrix::Diagonal(x = LD)
    Dinv_diag <- Matrix::Diagonal(x = 1/D)
    LDI_diag  <- Matrix::Diagonal(x = 1 + LD)
    
    # blocks
    H11 <- Matrix::t(B) %*% Dinv_diag %*% LD_diag %*% Dinv_diag %*% B
    H12 <- Matrix::t(B) %*% Dinv_diag %*% (dld_dw)
    H13 <- Matrix::t(B) %*% (LDI_diag %*% Dinv_diag) - Dinv_diag
    H14 <- Matrix::Matrix(0, nrow = n, ncol = 1, sparse = TRUE)
    H11 <- H11[-n, -n]
    H12 <- H12[-n, , drop = FALSE]
    H13 <- H13[-n, -n]
    H14 <- H14[-n, , drop = FALSE]
    first_row <- cbind(H11, H12, H13, H14)
    
    H21 <- Matrix::Matrix(Matrix::t(H12), sparse = TRUE)
    H22 <- Matrix::Matrix( -(nu - 1) / (w^2) +
                     sum((a/D)^2 * a_tilde_shift * pl -
                           (b/D)^2 * b_tilde_shift * pu -
                           (b/D * pu - a/D * pl)^2),
                   nrow = 1, ncol = 1, sparse = TRUE)
    H23 <- Matrix::Matrix(dld_dw, nrow = 1, sparse = TRUE)
    H23 <- H23[, -n, drop = FALSE]
    H24 <- Matrix::Matrix(-1, nrow = 1, ncol = 1, sparse = TRUE)
    second_row <- cbind(H21, H22, H23, H24)
    
    H31 <- Matrix::Matrix(Matrix::t(H13), sparse = TRUE)
    H32 <- Matrix::Matrix(dld_dw, nrow = n, ncol = 1, sparse = TRUE)
    H33 <- LDI_diag
    H34 <- Matrix::Matrix(0, nrow = n, ncol = 1, sparse = TRUE)
    H32 <- H32[-n, , drop = FALSE]
    H33 <- H33[-n, -n]
    H34 <- H34[-n, , drop = FALSE]
    third_row <- cbind(H31, H32, H33, H34)
    
    H41 <- Matrix::Matrix(0, nrow = 1, ncol = n-1, sparse = TRUE)
    H42 <- Matrix::Matrix(-1, nrow = 1, ncol = 1, sparse = TRUE)
    H43 <- Matrix::Matrix(0, nrow = 1, ncol = n-1, sparse = TRUE)
    H44 <- Matrix::Matrix(-kap * exp(stats::dnorm(kap, log = TRUE) -
                               TruncatedNormal::lnNpr(-kap, Inf)) -
                    exp(stats::dnorm(kap, log = TRUE) -
                          TruncatedNormal::lnNpr(-kap, Inf))^2 + 1,
                  nrow = 1, ncol = 1, sparse = TRUE)
    fourth_row <- cbind(H41, H42, H43, H44)
    
    J <- rbind(first_row, second_row, third_row, fourth_row)
  } else {
    J <- NULL
  }
  
  ## return
  if (deriv == "grad") return(grad)
  if (deriv == "jac")  return(J)
  list(grad = grad, jac = J)
}






