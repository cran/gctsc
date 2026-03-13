#' Marginal Models for Copula Time Series
#'
#' These functions define the marginal distributions used in copula-based count time series models.
#'
#' @name marginal.gctsc
#' @aliases
#' poisson.marg
#' negbin.marg
#' binom.marg
#' zip.marg
#' zib.marg
#' bbinom.marg
#' zibb.marg
#'
#' @description
#' The following marginal models are currently available:
#' \describe{
#'   \item{\code{poisson.marg(link = "log")}}{Poisson distribution.}
#'   \item{\code{negbin.marg(link = "log")}}{Negative binomial distribution.}
#'   \item{\code{binom.marg(link = "logit", size)}}{Binomial distribution with fixed number of trials.}
#'   \item{\code{bbinom.marg(link = "logit", size)}}{Beta-binomial with overdispersion.}
#'   \item{\code{zip.marg(link = "log")}}{Zero-inflated Poisson model.}
#'   \item{\code{zib.marg(link = "logit", size)}}{Zero-inflated binomial.}
#'   \item{\code{zibb.marg(link = "logit", size)}}{Zero-inflated beta-binomial with separate covariates for zero inflation.}
#' }
#'


#' @param link The link function for the mean (e.g., \code{"log"}, \code{"logit"}, \code{"identity"}).
#' @param size Number of trials (for binomial-type models).
#' @param lambda.lower Optional lower bounds on parameters.
#' @param lambda.upper Optional upper bounds on parameters.
#'
#' @details
#' Each marginal constructor returns an object of class \code{"marginal.gctsc"} which defines:
#' \itemize{
#'   \item \code{start}: a function to compute starting values.
#'   \item \code{npar}: number of parameters.
#'   \item \code{bounds}: truncation bounds on the latent Gaussian.
#' }
#'
#' These marginals are designed to work with \code{gctsc()} and its related methods.
#'
#' @return A list of class \code{"marginal.gctsc"} representing the marginal model.
#'
#' @seealso \code{\link{gctsc}}, \code{\link{predict.gctsc}}, \code{\link{arma.cormat}}
#'
#' @references
#' Cribari-Neto, F. and Zeileis, A. (2010). Beta regression in R. \emph{Journal of Statistical Software}, \strong{34}(2): 1–24.
#'
#' Ferrari, S.L.P. and Cribari-Neto, F. (2004). Beta regression for modeling rates and proportions. \emph{Journal of Applied Statistics}, \strong{31}(7): 799–815.
#'
#' Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics}, \strong{6}: 1517–1549.

#' @examples
#' poisson.marg(link = "identity")
#' zibb.marg(link = "logit", size = 24)
#'
#' @rdname marginal.gctsc
#'@export
#' @usage poisson.marg(link = "identity", lambda.lower = NULL, lambda.upper = NULL)
poisson.marg <- function(link = "identity", lambda.lower = NULL, lambda.upper = NULL) {


  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )


  obj <- list(
    start = function(y, x) {
      fit <- glm.fit(x = x, y = y, family = poisson(link))
      lambda <- fit$coefficients
      if (anyNA(lambda)) stop("NA detected in coefficient estimates. Check design matrix.")

      names(lambda) <- prefixed_names(x, "mu_")
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "poisson.marg()")
      lambda
    },

    npar = function(x) NCOL(x),

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      mu <- invlink(x %*% lambda)
      if (any(mu < 0)) stop("Negative mean detected. Use 'log' link or check predictors.")
      pdf <- dpois(y, mu)
      cdf <- ppois(y, mu)
      bounds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bounds$lower, bounds$upper)
    }

  )
  class(obj) <- "marginal.gctsc"
  obj
}


#' Binomial marginal model (supports y as vector or cbind(success, failure))
#'
#' @rdname marginal.gctsc
#'@export
#' @usage binom.marg(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL)
binom.marg <- function(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL) {
  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    "logit" = plogis,
                    stop("Unsupported link function")
  )



  obj <- list(
    start = function(y, x) {
      if (is.vector(y) || NCOL(y) == 1) {
        if (is.null(size)) stop("If y has 1 column, size must be provided.")
        y_bin <- cbind(y, size - y)
      } else {
        y_bin <- y
      }

      df <- data.frame(y = I(y_bin), x)
      fit <- glm(y ~ . - 1, family = binomial(link=link), data = df)
      lambda <- coef(fit)
      if (anyNA(lambda)) stop("NA in coefficient estimates. Check covariates for collinearity or data issues.")

      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "binom.marg()")

      names(lambda) <- prefixed_names(x, "mu_")
      lambda
    },

    npar = function(x) NCOL(x),

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      if (is.vector(y) || NCOL(y) == 1) {
        if (is.null(size)) stop("If y has 1 column, size must be provided.")
        size <- rep(size, length(y))
        successes <- y
      } else {
        size <- rowSums(y)
        successes <- y[, 1]
      }

      mu <- invlink(x %*% lambda)
      pdf <- dbinom(successes, size, mu)
      cdf <- pbinom(successes, size, mu)
      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}



#' Zero-Inflated Binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zib.marg(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL)
zib.marg <- function(link = "logit", size = NULL, lambda.lower = NULL, lambda.upper = NULL) {
  if (is.null(size)) stop("size must be provided for ZIB marginal.")

  invlink <- switch(link,
                    "logit" = plogis,
                    stop("Unsupported link function"))

  obj <- list(
    start = function(y, x) {
      X_mu  <- x$mu
      X_pi0 <- x$pi0
      if (!is.numeric(y)) stop("y must be numeric counts.")
      
      size_vec <- rep(size, length(y))
      
      # Fit binomial only on non-zero observations
      nonzero_idx <- (y > 0)
      if (any(nonzero_idx)) {
        fit_mu <- glm(cbind(y, size_vec - y) ~ X_mu - 1,
                      family = binomial(link = link),
                      subset = nonzero_idx)
        beta <- coef(fit_mu)
      } else {
        beta <- rep(0, ncol(X_mu))  # fallback
      }
      
      mu <- invlink(X_mu %*% beta)
      
      # Fit zero-inflation model (logit link)
      f0 <- dbinom(0, size_vec, mu)
      
      # Crude estimate of structural zero probability
      pi0_est <- pmax((as.numeric(y == 0) - f0) / (1 - f0), 1e-6)
      
      # Fit α on this adjusted probability
      fit_pi0 <- glm(pi0_est ~ X_pi0 - 1,
                     family = quasibinomial(link = "logit"))
      # fit_pi0 <- glm(I(y == 0) ~ X_pi0 - 1,
      #                family = binomial(link = "logit"))
      alpha <- coef(fit_pi0)
      
      lambda <- c(beta, alpha)
      names(lambda) <- c(prefixed_names(X_mu, "mu_"),
                         prefixed_names(X_pi0, "pi0_"))
      
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zib.marg()")
      lambda
      },

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu) + ncol(x$pi0)
    },

    bounds = function(y, x, lambda, family = "gaussian", df=NULL) {
      if (!is.list(x)) stop("x must be a list with 'mu' and 'pi0'")
      X_mu <- x$mu
      X_pi0 <- x$pi0
      p_mu <- ncol(X_mu)
      beta <- lambda[1:p_mu]
      alpha <- lambda[(p_mu + 1):length(lambda)]

      mu <- invlink(X_mu %*% beta)
      pi0 <- plogis(X_pi0 %*% alpha)
      size <- rep(size, length(y))
      f0 <- dbinom(0, size, mu)
      pmf <- dbinom(y, size, mu)
      cdf <- pbinom(y, size, mu)

      pdf <- ifelse(y == 0,
                    pi0 + (1 - pi0) * f0,
                    (1 - pi0) * pmf)
      cdf <- pi0 + (1 - pi0) * cdf

      bds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bds$lower, bds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}


#' Negative binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage negbin.marg(link = "identity", lambda.lower = NULL, lambda.upper = NULL)
negbin.marg <- function(link = "identity" ,lambda.lower = NULL, lambda.upper = NULL) {
  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )

  obj <- list(
    start = function(y, x) {
      eps <- sqrt(.Machine$double.eps)
      m <- glm.fit(x, y, family = poisson(link=link))
      mu <- fitted(m)
      kappa <- max(10 * eps, mean(((y - mu)^2 - mu) / mu^2))
      lambda <- c(coef(m), dispersion = kappa)
      xnames <- .default_colnames(x, "X")
      names(lambda) <- c(paste0("mu_", xnames), "dispersion")


      # Check bounds
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "negbin.marg()")
      lambda
    }
    ,

    npar = function(x) (NCOL(x) +1) ,

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      beta <- lambda[1:NCOL(x)]
      dispersion <- lambda[length(lambda)]
      mu <- invlink(x %*% beta)
      if (any(mu < 0)) stop("Negative mean under identity link. Consider using 'log' link.")
      size <- 1 / dispersion
      pdf <- dnbinom(y, mu = mu, size = size)
      cdf <- pnbinom(y, mu = mu, size = size)
      bounds <- safe_cdf_bounds(pdf, cdf,family = family, df = df)
      cbind(bounds$lower, bounds$upper)
    }
    

  )
  class(obj) <- "marginal.gctsc"
  obj
}

#' Zero-Inflated Poisson marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zip.marg(link = "identity", lambda.lower = NULL, lambda.upper = NULL)
zip.marg <- function(link = "identity", lambda.lower = NULL, lambda.upper = NULL) {
  invlink <- switch(link,
                    "identity" = function(eta) eta,
                    "log" = exp,
                    stop("Unsupported link function")
  )

  obj <- list(
    start = function(y, x) {

      X_mu <- x$mu
      X_pi0 <- x$pi0
      if (!is.numeric(y)) stop("y must be numeric.")

      fit <- glm.fit(X_mu, y, family = poisson(link = link))
      beta <- coef(fit)
      mu <- invlink(X_mu %*% beta)

      # Estimate zero-inflation
      f0 <- dpois(0, mu)
      pi0_hat <- mean(y == 0) - mean(f0)
      pi0_hat <- min(max(pi0_hat, 1e-6), 1 - 1e-6)

      fit_pi0 <- glm(I(y == 0) ~ X_pi0 - 1, family = binomial(link = "logit"))
      alpha <- coef(fit_pi0)

      lambda <- c(beta, alpha)
      names(lambda) <- c(prefixed_names(X_mu, "mu_"),
                         prefixed_names(X_pi0, "pi0_"))


      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zip.marg()")
      lambda
    },

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu) + ncol(x$pi0)
    },

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      X_mu <- x$mu
      X_pi0 <- x$pi0
      p_mu <- ncol(X_mu)
      beta <- lambda[1:p_mu]
      alpha <- lambda[(p_mu + 1):length(lambda)]

      mu <- invlink(X_mu %*% beta)
      pi0 <- plogis(X_pi0 %*% alpha)

      f0 <- dpois(0, mu)
      pmf <- dpois(y, mu)
      cdf <- ppois(y, mu)

      pdf <- ifelse(y == 0,
                    pi0 + (1 - pi0) * f0,
                    (1 - pi0) * pmf)
      cdf <- pi0 + (1 - pi0) * cdf

      bds <- safe_cdf_bounds(pdf, cdf, family, df)
      cbind(bds$lower, bds$upper)
      }
  )

  class(obj) <- "marginal.gctsc"
  obj
}


#' Beta-Binomial marginal model
#'
#' @rdname marginal.gctsc
#'@export
#' @usage bbinom.marg(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL)
bbinom.marg <- function(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL) {
  if (missing(size)) stop("size must be provided for beta-binomial marginal.")

  invlink <- switch(link,
                    "logit" = plogis,
                    stop("Unsupported link"))

  obj <- list(
    size = size,

    start = function(y, x) {
      assigns <- attr(x, "assign")
      has_intercept <- has_intercept(x)


      if (ncol(x) == 1 && has_intercept){
        df <- data.frame(y = y, trials = size)
        fit <- VGAM::vglm(cbind(y, trials - y) ~ 1, family = VGAM::betabinomial, data = df)
        cf <- VGAM::Coef(fit)
        beta <- qlogis(cf[1])
        rho_logit <- qlogis(cf[2])
        lambda <- c(beta, rho_logit)
        names(lambda) <- c("mu_(Intercept)", "logit_rho")

      } else {
        if (has_intercept) {
          x_sub <- x[, -1, drop = FALSE]
          coef_names <- c("mu_Intercept", prefixed_names(x_sub, "mu_"))
        } else {
          x_sub <- x
          coef_names <-   prefixed_names(x_sub, "mu_")
        }


        df <- data.frame(y = y, trials = size, x_sub)
        formula <- as.formula(paste("cbind(y, trials - y) ~", paste(colnames(x_sub), collapse = " + ")))
        fit <- VGAM::vglm(formula, family = VGAM::betabinomial, data = df)
        cf <- VGAM::Coef(fit)
        n_coef <- length(cf)
        beta <- c(cf[1], cf[3:(n_coef)])
        rho_logit <- (cf[2])

        lambda <- c(beta, rho_logit)
        
        names(lambda) <- c(coef_names, "logit_rho")
        lambda
        }



      # Set and validate bounds
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "bbinom.marg()")
      lambda
    },

    npar = function(x) NCOL(x) + 1,  # slopes + rho

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      beta <- lambda[1:NCOL(x)]
      rho_logit <- lambda[length(lambda)]
      rho <- plogis(rho_logit)
      mu   <- invlink(x %*% beta)

      alpha1 <- mu * (1 - rho) / rho
      alpha2 <- (1 - mu) * (1 - rho) / rho

      pdf <- VGAM::dbetabinom.ab(y, size = size, shape1 = alpha1, shape2 = alpha2, log = FALSE)
      cdf <- VGAM::pbetabinom.ab(y, size = size, shape1 = alpha1, shape2 = alpha2, log = FALSE)

      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}



#' Zero-Inflated Beta-Binomial with separate seasonal structure for pi0 and BB mean
#'
#' @rdname marginal.gctsc
#'@export
#' @usage zibb.marg(link = "logit", size,  lambda.lower = NULL, lambda.upper = NULL)
zibb.marg <- function(link = "logit", size, lambda.lower = NULL, lambda.upper = NULL) {
  if (missing(size)) stop("size must be provided.")

  invlink <- switch(link,
                    "logit" = plogis,
                    stop("Unsupported link function for mean")
  )
  obj <- list(
    size = size,

    start = function(y, x) {
      X_mu <- x$mu
      X_pi0 <- x$pi0
      has_intercept <- has_intercept(X_mu)
      y_nonzero <- y[y > 0]
      if (has_only_intercept(X_mu)) {
        df_nonzero <- data.frame(y = y_nonzero, trials = size)
        fit <- VGAM::vglm(cbind(y, trials - y) ~ 1,
                          family = VGAM::betabinomial,
                          data = df_nonzero)
        cf <- VGAM::Coef(fit)
        beta <- qlogis(cf[1])
        rho_logit <- qlogis(cf[2])
        coef_names <- c("mu_Intercept", "logit_rho")

      } else {
        if (has_intercept) {
          intercept_idx <- which(is_intercept_col(X_mu))
          X_mu_sub <- X_mu[, -intercept_idx, drop = FALSE]
          coef_names <- c("mu_Intercept",prefixed_names(X_mu_sub, "mu_"), "logit_rho")
        } else {
          X_mu_sub <- X_mu
          coef_names <-c( prefixed_names(X_mu_sub, "mu_"), "logit_rho")
        }


        df_nonzero <- data.frame(y = y_nonzero,
                                 trials = size,
                                 X_mu[y > 0, , drop = FALSE])
        formula <- as.formula(paste("cbind(y, trials - y) ~", paste(colnames(X_mu_sub), collapse = " + ")))
        fit <- VGAM::vglm(formula, family = VGAM::betabinomial, data = df_nonzero)
        cf <- VGAM::Coef(fit)
        n_coef <- length(cf)
        beta <- c(cf[1], cf[3:(n_coef)])
        rho_logit <-  (cf[2])
      }

      # Fit logistic regression for π₀
      fit_pi0 <- glm(I(y == 0) ~ X_pi0 - 1, family = binomial(link = "logit"))
      alpha <- coef(fit_pi0)
      lambda <- c(beta, rho_logit, alpha)
      names(lambda) <- c(coef_names, prefixed_names(X_pi0, "pi0_"))


      # Attach bounds if provided
      lambda <- attach_bounds(lambda, lambda.lower, lambda.upper, "zibb.marg()")
      lambda
    }

    ,

    npar = function(x) {
      if (!is.list(x) || is.null(x$mu) || is.null(x$pi0)) {
        stop("x must be a list with elements 'mu' and 'pi0'")
      }
      ncol(x$mu)+1 + ncol(x$pi0)
    },

    bounds = function(y, x, lambda,family = "gaussian", df=NULL) {
      if (!is.list(x)) stop("x must be a list with 'mu' and 'pi0'")
      X_mu <- x$mu
      X_pi0 <-x$pi0
      p_mu <- ncol(X_mu)
      np <- ncol(X_pi0)
      beta <- lambda[1:p_mu]
      rho_logit <-lambda[p_mu + 1]
      rho <- plogis(rho_logit)
      alpha <- lambda[(p_mu + 2):(p_mu + 1 + np)]


      mu <- invlink(X_mu %*% beta)
      pi0 <- plogis(X_pi0 %*% alpha)

      alpha1 <- mu * (1 - rho) / rho
      alpha2 <- (1 - mu) * (1 - rho) / rho

      f0 <- VGAM::dbetabinom.ab(0, size = size, shape1 = alpha1, shape2 = alpha2)
      pmf <- VGAM::dbetabinom.ab(y, size = size, shape1 = alpha1, shape2 = alpha2)
      cdf <- VGAM::pbetabinom.ab(y, size = size, shape1 = alpha1, shape2 = alpha2)

      pdf <- ifelse(y == 0, pi0 + (1 - pi0) * f0, (1 - pi0) * pmf)
      cdf <- pi0 + (1 - pi0) * cdf

      bounds <- safe_cdf_bounds(pdf, cdf,family, df)
      cbind(bounds$lower, bounds$upper)
    }
  )

  class(obj) <- "marginal.gctsc"
  obj
}

