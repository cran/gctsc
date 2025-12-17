
## -----------------------------------------------------------
## simulation Validate helpers 
## -----------------------------------------------------------

.scalar_or_nsim <- function(x, nsim, name, fn) {
  len <- length(x)
  if (len %in% c(1L, nsim)) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      stop(sprintf("%s(): '%s' must be numeric and finite.", fn, name), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  stop(sprintf("%s(): '%s' must be scalar or length nsim = %d. Got length %d.",
               fn, name, nsim, len), call. = FALSE)
}

.check_size <- function(size, fn, scalar_only = TRUE, nsim = NULL) {
  if (scalar_only) {
    if (length(size) != 1L)
      stop(sprintf("%s(): 'size' must be a single positive integer.", fn), call. = FALSE)
  } else {
    .scalar_or_nsim(size, nsim, "size", fn)
  }
  if (!is.numeric(size) || !is.finite(size) || size <= 0 || size != as.integer(size)) {
    stop(sprintf("%s(): 'size' must be a positive integer.", fn), call. = FALSE)
  }
}

.check_common <- function(nsim, tau, arma_order, seed, fn) {
  # nsim
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 1 || nsim != as.integer(nsim)) {
    stop(sprintf("%s(): 'nsim' must be a single positive integer. Got %s.",
                 fn, paste(nsim, collapse = ",")), call. = FALSE)
  }
  # arma_order
  if (!is.numeric(arma_order) || length(arma_order) != 2L ||
      any(!is.finite(arma_order)) || any(arma_order < 0) ||
      any(arma_order != as.integer(arma_order))) {
    stop(sprintf("%s(): 'arma_order' must be integer c(p, q) with p,q >= 0. Got %s.",
                 fn, paste(arma_order, collapse = ",")), call. = FALSE)
  }
  if (all(arma_order == 0L)) {
    stop(sprintf("%s(): ARMA(0,0) is not supported.", fn), call. = FALSE)
  }
  # tau
  if (!is.numeric(tau) || any(!is.finite(tau))) {
    stop(sprintf("%s(): 'tau' must be numeric and finite.", fn), call. = FALSE)
  }
  if (length(tau) != sum(arma_order)) {
    stop(sprintf("%s(): length(tau) must equal p+q = %d. Got %d.",
                 fn, sum(arma_order), length(tau)), call. = FALSE)
  }
  # seed
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed))) {
    stop(sprintf("%s(): 'seed' must be NULL or a single integer.", fn), call. = FALSE)
  }
  
  # Stationarity/invertibility checks
  p <- arma_order[1]; q <- arma_order[2]
  phi   <- if (p > 0) tau[seq_len(p)] else numeric(0)
  theta <- if (q > 0) tau[p + seq_len(q)] else numeric(0)
  
  if (p > 0) {
    ar_poly <- c(1, -phi)                     # 1 - sum phi_i z^i
    r <- polyroot(ar_poly)
    if (any(Mod(r) <= 1)) {
      stop(sprintf("%s(): AR polynomial is not stationary (a root has |z| <=1). Coefficients: %s",
                   fn, paste(phi, collapse = ",")), call. = FALSE)
    }
  }
  if (q > 0) {
    ma_poly <- c(1, theta)                     # 1 + sum theta_j z^j
    r <- polyroot(ma_poly)
    if (any(Mod(r) <= 1)) {
      stop(sprintf("%s(): MA polynomial is not invertible (a root has |z| <= 1). Coefficients: %s",
                   fn, paste(theta, collapse = ",")), call. = FALSE)
    }
  }
  invisible(TRUE)
}

.check_mu <- function(mu, nsim, fn) {
  .scalar_or_nsim(mu, nsim, "mu", fn)
  if (any(!(mu > 0))) {
    stop(sprintf("%s(): 'mu' must be strictly > 0.", fn), call. = FALSE)
  }
}

.check_prob <- function(prob, nsim, fn) {
  .scalar_or_nsim(prob, nsim, "prob", fn)
  if (any(!(prob > 0 & prob < 1))) {
    stop(sprintf("%s(): 'prob' must lie in (0, 1).", fn), call. = FALSE)
  }
}

.recyclen <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) == n) return(x)
  if (length(x) == 1L) return(rep(x, n))
  stop(sprintf("Parameter '%s' must be length 1 or %d (got %d).", name, n, length(x)))
}


.check_pi0 <- function(pi0, nsim, fn) {
  .scalar_or_nsim(pi0, nsim, "pi0", fn)
  if (any(!(pi0 >= 0 & pi0 < 1))) {
    stop(sprintf("%s(): 'pi0' must lie in [0, 1).", fn), call. = FALSE)
  }
}

.check_dispersion <- function(dispersion, nsim, fn) {
  .scalar_or_nsim(dispersion, nsim, "dispersion", fn)
  if (any(!(dispersion > 0))) {
    stop(sprintf("%s(): 'dispersion' must be strictly > 0.", fn), call. = FALSE)
  }
}


.check_rho <- function(rho, nsim, fn) {
  .scalar_or_nsim(rho, nsim, "rho", fn)
  if (any(!(rho > 0 & rho < 1))) {
    stop(sprintf("%s(): 'rho' must lie in (0, 1).", fn), call. = FALSE)
  }
}

# ----------------------------------
# Common validators for fitting
# -----------------------------------
# small helper
`%||%` <- function(a, b) if (is.null(a)) b else a

.validate_method <- function(method, fn) {
  ok <- c("TMET","GHK","CE","VMET")
  if (!method %in% ok) {
    stop(sprintf("%s(): 'method' must be one of %s. Got %s.",
                 fn, paste(ok, collapse = ", "), method), call. = FALSE)
  }
}

.validate_marg_cormat <- function(marginal, cormat, fn) {
  if (is.function(marginal)) marginal <- marginal()
  if (!inherits(marginal, "marginal.gctsc"))
    stop(sprintf("%s(): 'marginal' must be a 'marginal.gctsc' object.", fn), call. = FALSE)
  if (is.function(cormat)) cormat <- cormat()
  if (!inherits(cormat, "cormat.gctsc"))
    stop(sprintf("%s(): 'cormat' must be a 'cormat.gctsc' object (e.g., arma.cormat()).", fn), call. = FALSE)
  invisible(list(marginal = marginal, cormat = cormat))
}


.validate_options <- function(method, QMC, options, fn) {
  # options can be NULL
  if (is.null(options)) return(invisible(NULL))
  if (!is.list(options)) {
    stop(sprintf("%s(): 'options' must be a list or NULL.", fn), call. = FALSE)
  }
  
  # QMC must be logical scalar (for all methods)
  if (!is.logical(QMC) || length(QMC) != 1L || is.na(QMC)) {
    stop(sprintf("%s(): 'QMC' must be TRUE/FALSE.", fn), call. = FALSE)
  }
  
  # GHK / TMET: only validate M/seed if present
  if (!is.null(options$M)) {
    M <- options$M
    ok <- is.numeric(M) && length(M) == 1L && is.finite(M) && M == as.integer(M) && M > 0
    if (!ok) {
      stop(sprintf("%s(): options$M must be a single positive integer when supplied.", fn), call. = FALSE)
    }
  }
  if (!is.null(options$seed)) {
    s <- options$seed
    ok <- is.numeric(s) && length(s) == 1L && is.finite(s) && s == as.integer(s)
    if (!ok) {
      stop(sprintf("%s(): options$seed must be NULL or a single integer.", fn), call. = FALSE)
    }
  }
  
  invisible(NULL)
}


# For ZI vs non-ZI formula handling
.validate_formula_input <- function(formula, marginal, fn) {
  # Accept: a) a single formula, b) NULL/missing, c) a named list
  if (inherits(formula, "formula")) {
    formula <- list(mu = formula)
  } else if (missing(formula) || is.null(formula)) {
    formula <- list()
  } else if (!is.list(formula)) {
    stop(sprintf("%s(): 'formula' must be a formula or a named list of formulas.", fn), call. = FALSE)
  }
  
  # alias handling: allow 'pi' as 'pi0'
  if (!is.null(formula$pi) && is.null(formula$pi0)) {
    message(sprintf("%s(): 'pi' detected; treating as 'pi0'.", fn))
    formula$pi0 <- formula$pi
    formula$pi  <- NULL
  }
  
  if (has_ZI(marginal)) {
    # ZI needs mu and pi0; default to intercept-only if missing
    if (is.null(formula$mu) && is.null(formula$pi0)) {
      message(sprintf("%s(): zero-inflated marginal with no formulas; using mu ~ 1 and pi0 ~ 1.", fn))
      formula$mu  <- y ~ 1
      formula$pi0 <- ~ 1
    } else {
      if (is.null(formula$mu))  stop(sprintf("%s(): supply formula$mu (e.g., y ~ 1).", fn), call. = FALSE)
      if (is.null(formula$pi0)) {
        message(sprintf("%s(): 'pi0' not supplied; using pi0 ~ 1.", fn))
        formula$pi0 <- ~ 1
      }
    }
    if (!inherits(formula$mu,  "formula")) stop(sprintf("%s(): formula$mu must be a formula.",  fn), call. = FALSE)
    if (!inherits(formula$pi0, "formula")) stop(sprintf("%s(): formula$pi0 must be a formula.", fn), call. = FALSE)
    
  } else {
    # non-ZI: only mu
    if (is.null(formula$mu)) {
      if (length(formula) == 1L && inherits(formula[[1L]], "formula"))
        formula$mu <- formula[[1L]]
      else
        stop(sprintf("%s(): provide formula for non-ZI marginals (e.g., y ~ 1).", fn), call. = FALSE)
    }
    if (!inherits(formula$mu, "formula")) stop(sprintf("%s(): must be a formula.", fn), call. = FALSE)
    # ignore any pi/pi0 politely
    if (!is.null(formula$pi0) || !is.null(formula$pi)) {
      message(sprintf("%s(): ignoring pi/pi0 formula for non-zero-inflated marginal.", fn))
      formula$pi0 <- NULL; formula$pi <- NULL
    }
  }
  
  formula
}


.build_design <- function(formula, data, marginal, fn) {
  # Validate and normalize the formula parts
  formula <- .validate_formula_input(formula, marginal, fn)
  
  # Pick the mean formula
  f_mu <- formula$mu
  
  # === Resolve data if not provided ===
  if (missing(data) || is.null(data)) {
    vars <- all.vars(f_mu)
    y_name <- as.character(f_mu[[2L]])
    intercept_only <- length(attr(terms(f_mu), "term.labels")) == 0L
    
    # Helper to fetch response from env
    get_response <- function(name) {
      get0(name, envir = environment(f_mu), inherits = TRUE) %||%
        get0(name, envir = parent.frame(), inherits = TRUE)
    }
    
    y_val <- get_response(y_name)
    if (is.null(y_val)) {
      stop(sprintf("%s(): could not find response '%s'. Pass data=... or bind it in the calling env.",
                   fn, y_name), call. = FALSE)
    }
    if (!intercept_only && length(vars) > 1L) {
      stop(sprintf("%s(): please supply 'data=...' when using covariates.", fn), call. = FALSE)
    }
    data <- setNames(data.frame(y_val), y_name)
  }
  
  # === Build model frames ===
  if (has_ZI(marginal)) {
    mf_mu  <- model.frame(f_mu,            data = data, na.action = na.pass)
    mf_pi0 <- model.frame(formula$pi0,     data = data, na.action = na.pass)
    y      <- model.response(mf_mu)
    if (!is.numeric(y)) {
      stop(sprintf("%s(): response must be numeric counts.", fn), call. = FALSE)
    }
    X_mu   <- model.matrix(f_mu,  mf_mu)
    X_pi0  <- model.matrix(formula$pi0, mf_pi0)
    return(list(
      y = y,
      x = list(mu = X_mu, pi0 = X_pi0),
      terms = terms(f_mu),
      model = mf_mu
    ))
  } else {
    mf <- model.frame(f_mu, data = data, na.action = na.pass)
    y  <- model.response(mf)
    if (!is.numeric(y)) {
      stop(sprintf("%s(): response must be numeric counts.", fn), call. = FALSE)
    }
    X  <- model.matrix(f_mu, mf)
    return(list(
      y = y,
      x = X,
      terms = terms(f_mu),
      model = mf
    ))
  }
}


# x structure and nrow checks used in both gctsc() and gctsc.fit()
validate_x_structure <- function(x, marginal, fn = "gctsc") {
  if (has_ZI(marginal)) {
    if (!is.list(x) || is.null(x$mu))
      stop(sprintf("%s(): x must be list(mu = X_mu, pi0 = X_pi0) for zero-inflated marginals.", fn), call. = FALSE)
    if (!is.matrix(x$mu))  x$mu  <- as.matrix(x$mu)
    if (is.null(x$pi0)) x$pi0 <- matrix(1, nrow = nrow(x$mu), ncol = 1L)
    if (!is.matrix(x$pi0)) x$pi0 <- as.matrix(x$pi0)
  } else {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.matrix(x)) stop(sprintf("%s(): x must be a matrix for non-zero-inflated marginals.", fn), call. = FALSE)
  }
  invisible(TRUE)
}

check_x_nrow_matches_y <- function(x, y, marginal = NULL, fn = "gctsc") {
  n <- length(y)
  if (is.list(x)) {
    nm <- vapply(x, nrow, integer(1))
    if (any(nm != n))
      stop(sprintf("%s(): nrow mismatch: y has %d rows; got nrow(mu)=%d, nrow(pi0)=%d.",
                   fn, n, nrow(x$mu), nrow(x$pi0)), call. = FALSE)
  } else {
    if (nrow(x) != n)
      stop(sprintf("%s(): nrow(x) must equal length(y). Got %d vs %d.",
                   fn, nrow(x), n), call. = FALSE)
  }
  invisible(TRUE)
}

# Optional AR/MA admissibility check for initial tau
.check_arima_admissibility <- function(tau, p, q, fn) {
  if (p + q == 0) return(invisible(TRUE))
  if (!is.numeric(tau) || any(!is.finite(tau)) || length(tau) != (p + q))
    stop(sprintf("%s(): initial tau must be numeric, finite, length p+q.", fn), call. = FALSE)
  phi   <- if (p > 0) tau[seq_len(p)] else numeric(0)
  theta <- if (q > 0) tau[p + seq_len(q)] else numeric(0)
  if (p > 0) {
    ar_poly <- c(1, -phi)
    if (any(Mod(polyroot(ar_poly)) <= 1))
      stop(sprintf("%s(): AR polynomial not stationary at initial tau.", fn), call. = FALSE)
  }
  if (q > 0) {
    ma_poly <- c(1, theta)
    if (any(Mod(polyroot(ma_poly)) <= 1))
      stop(sprintf("%s(): MA polynomial not invertible at initial tau.", fn), call. = FALSE)
  }
  invisible(TRUE)
}

# --- Ensure intercept column exists if required ---
add_intercept_if_needed <- function(X_test, req_cols) {
  if ("(Intercept)" %in% req_cols && !"(Intercept)" %in% names(X_test)) {
    X_test <- c(`(Intercept)` = 1, X_test)
  }
  # Ensure it is a named numeric vector
  X_test <- unlist(X_test, use.names = TRUE)
  return(X_test)
}

