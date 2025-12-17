#' @useDynLib gctsc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
# core_utils.R
# Internal utilities for TMET and GHK simulation, gradient, and likelihood computation

#' @keywords internal
#' @noRd
sample_latent <- function(obj, lower, upper, delta = NULL, M = 1000, QMC = TRUE, ret_llk = TRUE,
                          method = c("TMET", "GHK")) {
  method <- match.arg(method)

  sampler_input <- list(
    a = lower,
    b = upper,
    condSd = sqrt(obj$cond_var),
    M = M,
    phi = obj$phi,
    q = obj$q,
    m = obj$m,
    Theta = obj$Theta,
    QMC = QMC
  )

  # Include delta only for TMET


  if (ret_llk) {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = ptmvn_tmet(sampler_input),
           "GHK"  = ptmvn_ghk(sampler_input),
           stop("Unknown method")
    )
  } else {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = rtmvn_tmet(sampler_input),
           "GHK"  = rtmvn_ghk(sampler_input),
           stop("Unknown method")
    )
  }
}

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

