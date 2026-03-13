#' Extract Coefficients from a gctsc Model
#'
#' Returns the estimated coefficients from a fitted `gctsc` model object.
#'
#' @param object An object of class `gctsc`.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @return A named numeric vector of model coefficients.
#' @export
coef.gctsc <- function(object,...) object$coef



#' Print a gctsc Model Object
#'
#' Displays the call, estimation method, parameter estimates, and likelihood information.
#'
#' @param x An object of class `gctsc`.
#' @param digits Number of significant digits to display.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @export
print.gctsc <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n--- Copula Count Time Series Model ---\n")
  cat("Copula family:", x$family, "\n")
  if (!is.null(x$df)) {
    cat("Degrees of freedom:", x$df, "\n")
  }
  cat("Sample size (n):", x$n, "\n")

  method_label <- switch(x$method,
                         "TMET" = "TMET: Tilted Importance Sampling",
                         "GHK"  = "GHK: Sequential Importance Sampling",
                         "CE"   = "CE: Continuous Extension Approximation",
                         paste("Unknown method:", x$method))
  cat("Estimation method:", method_label, "\n\n")

  cat("Call:\n", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "\n\n")

  if (x$convergence != 0) {
    cat("Warning: model did not converge (code", x$convergence, ")\n\n")
  }

  coefs <- x$coef
  se <- x$se
  has_se <- !is.null(se) && length(se) == length(coefs)

  if (x$nbeta > 0) {
    cat("Marginal parameters (psi):\n")
    beta_vals <- coefs[x$ibeta]
    if (has_se) {
      beta_se <- se[x$ibeta]
      printCoefmat(cbind(Estimate = beta_vals, `Std. Error` = beta_se), digits = digits)
    } else {
      print.default(format(beta_vals, digits = digits), print.gap = 2, quote = FALSE)
    }
    cat("\n")
  } else {
    cat("No marginal parameters estimated.\n\n")
  }

  if (x$ntau > 0) {
    cat("Copula dependence parameters (tau):\n")
    tau_vals <- coefs[x$itau]
    if (has_se) {
      tau_se <- se[x$itau]
      printCoefmat(cbind(Estimate = tau_vals, `Std. Error` = tau_se), digits = digits)
    } else {
      print.default(format(tau_vals, digits = digits), print.gap = 2, quote = FALSE)
    }
    cat("\n")
  } else {
    cat("No copula dependence parameters estimated.\n\n")
  }

  cat("Log-likelihood (approximate):", formatC(-x$maximum, digits = digits, format = "f"), "\n")
  cat("Monte Carlo draws (M):", x$options$M, "\n")

  invisible(x)
}

#' Summarize a gctsc Model Fit
#'
#' Computes standard errors, z-values, and p-values for the estimated parameters
#' in a fitted `gctsc` object.
#'
#' @param object An object of class `gctsc`.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @return A list of class `summary.gctsc` containing model summary statistics.
#' @export
summary.gctsc <- function(object, ...) {
  se <- object$se
  coef <- object$coef
  nbeta <- object$nbeta
  ntau <- object$ntau
  
  has_se <- !is.null(se) && is.numeric(se) && length(se) == length(coef)
  
  marginal_coef <- copula_coef <- NULL
  
  if (nbeta > 0) {
    beta <- coef[object$ibeta]
    if (has_se) {
      beta_se <- se[object$ibeta]
      zval <- beta / beta_se
      pval <- 2 * pnorm(-abs(zval))
    } else {
      beta_se <- zval <- pval <- rep(NA_real_, length(beta))
    }
    marginal_coef <- cbind(
      Estimate = beta,
      `Std. Error` = beta_se,
      `z value` = zval,
      `Pr(>|z|)` = pval
    )
    rownames(marginal_coef) <- names(beta)
  }
  
  if (ntau > 0) {
    tau <- coef[object$itau]
    if (has_se) {
      tau_se <- se[object$itau]
      zval <- tau / tau_se
      pval <- 2 * pnorm(-abs(zval))
    } else {
      tau_se <- zval <- pval <- rep(NA_real_, length(tau))
    }
    copula_coef <- cbind(
      Estimate = tau,
      `Std. Error` = tau_se,
      `z value` = zval,
      `Pr(>|z|)` = pval
    )
    rownames(copula_coef) <- names(tau)
  }
  
  loglik <- -object$maximum
  k <- length(coef)
  n <- object$n
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k
  
  structure(list(
    call = object$call,
    convergence = object$convergence,
    coefficients = list(marginal = marginal_coef, copula = copula_coef),
    loglik = loglik,
    aic = aic,
    bic = bic,
    maximum = -loglik
  ), class = "summary.gctsc")
}

#' Print Summary of a gctsc Model
#'
#' Displays summary statistics and model fit information for a fitted `gctsc` model.
#'
#' @param x An object of class `summary.gctsc`.
#' @param digits Number of significant digits to display.
#' @param ... Ignored. Included for S3 method compatibility.
#'
#' @export
print.summary.gctsc <- function(x, digits = 4, ...) {
  cat("\n--- Summary of Copula Time Series Model ---\n")
  cat("Call:\n", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "\n")

  if (x$convergence != 0) {
    cat("\nWarning: model did not converge (code", x$convergence, ")\n\n")
  }

  if (!is.null(x$coefficients$marginal)) {
    cat("\nMarginal Model Coefficients:\n")
    printCoefmat(x$coefficients$marginal, digits = digits, signif.legend = FALSE)
  } else {
    cat("\nNo coefficients in the marginal model.\n")
  }

  if (!is.null(x$coefficients$copula)) {
    cat("\nCopula (Dependence) Coefficients:\n")
    printCoefmat(x$coefficients$copula, digits = digits, signif.legend = FALSE)
  } else {
    cat("\nNo coefficients in the copula dependence model.\n")
  }

  # Safe signif. codes check
  all_coefs <- do.call(rbind, x$coefficients)
  if (!is.null(all_coefs) && ncol(all_coefs) >= 4) {
    pvals <- all_coefs[, 4L]
    if (getOption("show.signif.stars") &&
        any(!is.na(pvals) & pvals < 0.1)) {
      cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    }
  }

  cat("\nModel Fit Statistics:\n")
  cat("  Log-likelihood:", formatC(x$loglik, format = "f", digits = digits), "\n")
  cat("  AIC:            ", formatC(x$aic, format = "f", digits = digits), "\n")
  cat("  BIC:            ", formatC(x$bic, format = "f", digits = digits), "\n")

  invisible(x)
}
