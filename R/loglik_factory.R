#' @keywords internal
#' @noRd
build_loglik <- function(md, M, penalt=NA){
  n <- md$n
  y <- md$y
  x <- md$x
  is.int <- md$method !="CE"
  bounds <- md$marginal$bounds
  ibeta <- md$ibeta
  itau<- md$itau
  od <- md$cormat$od
  p <- md$cormat$od[1]  # AR order
  q <- md$cormat$od[2]  # MA order
  ns <- md$options$M
  seed <- md$options$seed
  c <- md$c
  method <- md$method
  QMC <- md$QMC
  pm <- md$pm
  cfg <- list(
    method = method,
    arg2   = if (is.int) ns else c,
    ret_llk = TRUE,
    pm = pm,
    od=od,
    QMC=QMC
  )
  eta <- md$fixed
  cache <- new.env()
  function( eta ) {
    beta <- eta[ibeta]
    if (!identical(cache$beta,beta)) {
      ab <-  bounds(y,x,beta)
      assign("beta",beta,envir=cache)
      assign("ab",ab,envir=cache)
    } else {
      ab <- get("ab",envir=cache)
    }
    if (is.null(ab) || any(is.nan(ab))) {


      return(NA)

    }
    tau <- eta[itau]
    if (!identical(cache$tau,tau)) {
      if (p > 0) {
        ar_coefs <- tau[1:p]  # First p elements are AR coefficients
        ar_roots <- polyroot(c(1, -ar_coefs))  # Note the negation for AR polynomial
        if (any(Mod(ar_roots) <= 1.01)) return(NA)  # Penalize invalid AR
      }

      # MA roots check (if q > 0)
      if (q > 0) {
        ma_coefs <- tau[(p + 1):(p + q)]  # Next q elements are MA coefficients
        ma_roots <- polyroot(c(1, ma_coefs))
        if (any(Mod(ma_roots) <= 1.01)) return(NA)  # Penalize invalid MA
      }
      assign("tau",tau,envir=cache)

    }
    if (is.int) set.seed(seed)
    lk <-  llk.fn(cfg, ab, tau)
    if ( is.finite(lk))   (lk) else penalt
  }
}


#' @keywords internal
#' @noRd
llk.fn <- function(cfg, ab, tau) {
  method <- cfg$method
  arg2 <- cfg$arg2
  ret_llk <- cfg$ret_llk
  od <- cfg$od
  QMC <- cfg$QMC
  pm <- cfg$pm
  result <- switch(method,
                   "CE"   = loglik_ce(ab, tau, c = arg2, od, ret_llk),
                   "GHK"  = loglik_ghk(ab, tau, M = arg2, od, QMC,ret_llk),
                   "TMET" = loglik_tmet(ab, tau, M = arg2, od,pm, QMC, ret_llk),
                   "VMET" = loglik_vmet(ab, tau,  M = arg2, od, QMC),
                   stop("Unknown method")
  )

  return(result)
}







