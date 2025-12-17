#' @keywords internal
#' @noRd
loglik_vmet <- function(ab, tau, M = 1000, od, QMC = TRUE, ret_llk = TRUE) {
  if (length(tau) != sum(od)) {
    stop("Length of 'tau' must equal sum of AR and MA orders: length(tau) = ",
         length(tau), ", expected = ", sum(od), ".")
  }
  if (all(od == 0)) {
    stop("ARMA(0,0) (white noise) is not supported. Please specify at least one AR or MA term.")
  }

  if (!is.matrix(ab) || ncol(ab) != 2) stop("ab must be an n x 2 matrix of bounds.")
  if (any(is.na(ab)) || any(ab[,1] > ab[,2])) {
    warning("Invalid bounds provided to VMET.")
    return(-1e20)
  }

  p <- od[1]
  q <- od[2]
  iar <- if (p) 1:p else NULL
  ima <- if (q) (p+1):(p+q) else NULL
  phi <- tau[iar]
  theta <- tau[ima]
  n <- nrow(ab)

  covM <- toeplitz(stats::ARMAacf(ar = phi, ma = theta, lag.max = n - 1, pacf = FALSE))

  result <- tryCatch({
    VeccTMVN::pmvn(lower = ab[,1], upper = ab[,2], sigma = covM,
                   reorder = 0, NLevel1 = 1, mean = rep(0, n),
                   NLevel2 = 1000, verbose = FALSE, retlog = TRUE, m = 30)
  }, error = function(e) {
    message("VMET failed: ", e$message)
    return(-1e20)
  })


  return(result)
}
