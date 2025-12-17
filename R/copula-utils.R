#' Set Options for Gaussian Copula Time Series Model
#'
#' Creates a control list for simulation and likelihood approximation in the
#' Gaussian copula model, including the random seed and Monte Carlo settings.
#'
#' @param seed Integer. Random seed for reproducibility (default: a random integer between 1 and 100000).
#' @param M Integer. Number of Monte Carlo samples used in the likelihood approximation (default: 1000).
#' @param ... Ignored. Included for S3 method compatibility.
#'
#'
#' @return
#' A list with components:
#'   \item{\code{seed}}{ Integer. The random seed used.}
#'   \item{\code{M}}{ Integer. Number of Monte Carlo samples.}
#'   \item{\code{opt}}{ A function used internally by \code{gctsc()} to
#'         perform optimization of the approximate log-likelihood.}

#' 
#' @export
gctsc.opts <- function(seed=NULL, M=1000, ...) {


  control <- list(...)
  opt <- function(start, llk_fn, lower, upper) {
    fn.opt <- function(x)  {if( any(x <= lower | x >= upper) ) (1e10)
      # Compute the log-likelihood sum
      # else {-llk_fn(x)}
      else {
        val <- -llk_fn(x)
        if (!is.finite(val)) return(1e10)
        return(val)
      }
    }
    ans <- optim(start, fn.opt, method= "BFGS", hessian = TRUE)
    if(ans$convergence) warning(paste("optim exits with code",ans$convergence))
    list(
      estimate = ans$par,
      maximum = ans$value,
      convergence = ans$convergence,
      hessian = ans$hessian
    )
  }
  list(seed=seed,M=M,opt=opt)
}


#' @keywords internal
#' @noRd
validate_x_structure <- function(x, marginal) {
  is_ZI <- has_ZI(marginal)

  if (is_ZI) {
    if (!is.list(x) || is.data.frame(x)) {
      stop("For zero-inflated marginals, `x` must be a list with elements `mu` and `pi0`.")
    }
    if (!all(c("mu", "pi0") %in% names(x))) {
      stop("For zero-inflated marginals, `x` must contain both `mu` and `pi0`.")
    }
  } else {
    if (is.list(x) && !is.data.frame(x)) {
      stop("For non-zero-inflated marginals, `x` must be a matrix.")
    }
  }
}



#' @keywords internal
#' @noRd
check_x_nrow_matches_y <- function(x, y) {
  if (!is.null(x)) {
    if (is.list(x)) {
      for (name in names(x)) {
        if (NROW(x[[name]]) != length(y)) {
          stop(sprintf("Design matrix '%s' must have the same number of rows as y", name))
        }
      }
    } else {
      if (NROW(x) != length(y)) {
        stop("Design matrix x must have the same number of rows as y")
      }
    }
  }
}


#' @keywords internal
#' @noRd
remove_na <- function(y, x) {
  if (is.null(x)) {
    x <- rep(1, length(y))  # intercept-only model
  }

  x_mat <- if (is.list(x)) do.call(cbind, x) else x
  not.na <- rowSums(is.na(cbind(y, x_mat))) == 0

  y_clean <- y[not.na]
  if (is.list(x)) {
    x_clean <- lapply(x, function(col) col[not.na])
  } else {
    x_clean <- x[not.na, , drop = FALSE]
  }

  list(y = y_clean, x = x_clean, not.na = not.na)
}
