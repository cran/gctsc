cg_solve <- function(Hfun, b, Mdiag = NULL, tol = 1e-6, maxiter = 1000, verbose = FALSE) {
  n <- length(b)
  x <- rep(0, n)        # initial guess
  r <- b - Hfun(x)      # residual
  
  # Jacobi preconditioner: diagonal of A
  if (is.null(Mdiag)) {
    z <- r              # no preconditioning
  } else {
    z <- r / Mdiag
  }
  
  p <- z
  rzold <- sum(r * z)
  
  for (i in seq_len(maxiter)) {
    Ap <- Hfun(p)
    alpha <- rzold / sum(p * Ap)
    
    x <- x + alpha * p
    r <- r - alpha * Ap
    
    if (sqrt(sum(r * r)) < tol) {
      if (verbose) cat(sprintf("   CG converged at iter %d: |r|=%.3e\n", i, sqrt(sum(r*r))))
      break
    }
    
    if (is.null(Mdiag)) {
      z <- r
    } else {
      z <- r / Mdiag
    }
    
    rznew <- sum(r * z)
    beta <- rznew / rzold
    p <- z + beta * p
    rzold <- rznew
    
    if (verbose) {
      cat(sprintf("   CG iter %d: |r|=%.3e\n", i, sqrt(sum(r*r))))
    }
  }
  x
}


lm_sparse_solver <- function(x0, jac, Condmv_Obj, a, b, nu,
                             maxit=50, tol=1e-8,
                             lambda0=1e-3, verbose=FALSE) {
  x <- x0
  lambda <- lambda0
  
  
  check_feasible <- function(z_delta) {
    n <- length(a)
    z <- delta <- rep(0, n)
    z[-n] <- z_delta[1:(n-1)]
    w <- z_delta[n]
    delta[-n] <- z_delta[(n+1):(2*n-1)]
    kap <- z_delta[2*n]
    
    a_scaled <- a / sqrt(nu)
    b_scaled <- b / sqrt(nu)
    B <- Condmv_Obj$B
    D <- sqrt(Condmv_Obj$cond_var)
    mu_c <- as.vector(B %*% z)
    
    a_tilde_shift <- (a_scaled * w - mu_c) / D - delta
    b_tilde_shift <- (b_scaled * w - mu_c) / D - delta
    
    bad <- which(a_tilde_shift >= b_tilde_shift - 1e-6 * abs(D) | w < 0)
    list(ok = length(bad) == 0, bad = bad)
  }
  
  
  
  
  for (k in seq_len(maxit)) {
    gj <- grad_jac_psiT(x, Condmv_Obj, a, b, nu, deriv="both")
    f  <- gj$grad      # gradient vector (∇Φ = Jᵀψ)
    J  <- gj$jac       # Jacobian (sparse)
    
    
    
    # g <- as.vector(crossprod(J, f))
    g <- as.vector(Matrix::t(J) %*% f)
    
    if (sqrt(sum(g^2)) < tol) {
      if (verbose) cat("Converged at iter", k, "\n")
      return(list(x=x, fval=f, iter=k, lambda=lambda, converged=TRUE))
    }
    
    # define Hfun for CG: (J^T J + lambda I) v
    Hfun <- function(v) { as.numeric(Matrix::t(J) %*% (J %*% v) + lambda * v) }
    #
    Mdiag <- Matrix::colSums(J^2) + lambda
    
    rhs <- -g
    step <- cg_solve(Hfun, rhs, Mdiag = Mdiag, tol = 1e-6, maxiter = 500, verbose = FALSE)
    
    # backtracking line search with feasibility
    t <- 1
    repeat {
      chk <- check_feasible(x + t*step)
      if (chk$ok) break
      if (verbose) cat("Infeasible indices at iter", k, ":", chk$bad, "\n")
      t <- t/2
      if (t < 1e-8) break
    }
    
    x_new <- x + t*step
    f_new <- grad_jac_psiT(x_new, Condmv_Obj, a, b, nu, deriv="grad")
    rho <- sum(f^2) - sum(f_new^2)
    
    if (rho > 0) {
      x <- x_new
      lambda <- lambda / 2
    } else {
      lambda <- lambda * 2
    }
    
    if (verbose) {
      cat(sprintf("Iter %d: |grad|=%.3e, step=%.2e, lambda=%.3e\n",
                  k, sqrt(sum(g^2)), norm(step, "2"), lambda))
    }
  }
  
  list(x=x, iter=maxit, lambda=lambda, converged=FALSE)
}