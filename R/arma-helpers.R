#' @keywords internal
aacvf <- function(a,h) {
  phi = a$phi
  theta = a$theta
  sigma2 = (1/(sum(ma.inf(a)^2)))
  p = length(phi)
  q = length(theta)

  m <- max(p, q + 1)
  gamma = numeric(max(p+1,h+1))

  psi = ma.inf(a,q)
  theta = c(1,theta)
  # For 0 <= k < m

  r = numeric(max(p+1,q+1,h+1))
  for (k in 0:(m - 1)) {
    limit <- max(0, q - k)
    r[k+1] <- sigma2*sum(theta[(k + 1):(k + limit + 1)] * psi[1:(limit + 1)])
  }

  first_col <- c(1, -phi)
  A <- toeplitz(first_col)
  gamma[1:(p+1)] = solve.qr(qr(A),r[1:(p+1)])

  #For k >= m, use pure AR recursion
  for (k in (p+1):h) {
    gamma[k + 1] <-  r[k+1] + sum(phi * gamma[k :(k - p +1 )])
  }

  return(gamma)
}

#' @keywords internal
ma.inf <- function(a, n = 200) {
  if (n == 0) return(1)
  theta <- c(a$theta, numeric(n - length(a$theta)))
  phi <- a$phi
  p <- length(phi)
  psi <- c(numeric(p), 1, numeric(n))
  idx_temp <- p:1
  for (j in 1:n) {
    idx <- idx_temp + j
    psi[j + p + 1] <- theta[j] + sum(phi * psi[idx])
  }
  return(psi[(0:n) + p + 1])
}

#' @keywords internal
ar.inf = function(a,n=100) {
  if (n == 0)
    return(1)
  phi = c(a$phi,numeric(n))
  theta = a$theta
  q = length(theta)
  pie = c(numeric(q),1,numeric(n))
  for (j in 1:n)
    pie[j+q+1] = -phi[j] - sum(theta * pie[(q:1)+j])
  return(pie[(0:n)+q+1])
}

#' @keywords internal
durbin_levinson <- function(gamma, n,v) {
  phi <- matrix(0, n, n)


  phi[1,1] <- gamma[2] / gamma[1]

  for (k in 2:n) {
    num <- gamma[k+1] - sum(phi[k-1, 1:(k-1)] * gamma[k:2])
    phi[k,k] <- num / v[k-1]

    for (j in 1:(k-1)) {
      phi[k,j] <- phi[k-1,j] - phi[k,k] * phi[k-1,k-j]
    }

  }

  # Return last row of phi as BLUP coefficients
  return(phi)
}

#' @keywords internal
check_tau_length <- function(tau, arma_order) {
  expected <- sum(arma_order)
  if (length(tau) != expected) {
    stop(sprintf("Invalid 'tau': expected length %d but got %d.", expected, length(tau)))
  }
}

