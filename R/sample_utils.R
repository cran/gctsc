#' @useDynLib gctsc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
# core_utils.R
# Internal utilities for TMET and GHK simulation, gradient, and likelihood computation

#' @keywords internal
#' @noRd
sample_mvn <- function(obj, lower, upper, delta = NULL, M = 1000, QMC = TRUE, ret_llk = TRUE,
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



#' @keywords internal
#' @noRd
sample_mvt <- function(obj, lower, upper, delta = NULL, M = 1000, QMC = TRUE, ret_llk = TRUE,
                       method = c("TMET", "GHK","GHK_MVT"), df) {
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
    QMC = QMC,
    df= df
    
  )
  
  # Include delta only for TMET
  
  
  if (ret_llk) {
    if (method == "TMET") {
      sampler_input$delta <- delta
    }
    switch(method,
           "TMET" = ptmvt_tmet(sampler_input),
           "GHK"  = ptmvmn_ghk(sampler_input),
           "GHK_MVT" = ptmvt_ghk(sampler_input),
           stop("Unknown method")
    )
  } else {
    rtmvt(sampler_input)
    
  }
  
  
}




