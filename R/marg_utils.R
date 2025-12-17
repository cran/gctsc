#' @keywords internal
#' @noRd
safe_cdf_bounds <- function(pdf, cdf) {
  EPS <- sqrt(.Machine$double.eps)
  EPS1 <- 1 - EPS
  cdf_l <- pmax(EPS, pmin(EPS1, cdf - pdf))
  cdf_u <- pmax(EPS, pmin(EPS1, cdf))
  cdf_l <- pmin(cdf_l, cdf_u - EPS)
  list(lower = qnorm(cdf_l), upper = qnorm(cdf_u))
}

#' @keywords internal
#' @noRd
has_ZI <- function(marginal) {
  fn_body <- body(marginal$start)
  any(grepl("x\\$pi", deparse(fn_body)))
}

#' @keywords internal
#' @noRd
has_intercept <- function(x) {
  assigns <- attr(x, "assign")

  if (!is.null(assigns)) {
    return(0 %in% assigns)
  }

  # Fallbacks if `assign` attribute is missing
  has_named_intercept <- !is.null(colnames(x)) &&
    any(tolower(colnames(x)) %in% c("intercept", "(intercept)"))

  has_column_of_ones <- any(apply(x, 2, function(col) all(col == 1)))

  return(has_named_intercept || has_column_of_ones)
}

#' @keywords internal
#' @noRd
has_only_intercept <- function(mat) {
  is.matrix(mat) && ncol(mat) == 1 && all(unique(mat) == 1)
}



is_intercept_col <- function(x) {
  apply(x, 2, function(col) all(col == 1))
}


# If X has no colnames, make X1, X2, ...
.default_colnames <- function(X, base = "X") {
  cn <- colnames(X)
  if (is.null(cn)) paste0(base, seq_len(ncol(X))) else cn
}

# Prefix names with a part label: "mu_" or "pi0_"
prefixed_names <- function(X, prefix) {
  paste0(prefix, .default_colnames(X, base = "X"))
}


attach_bounds <- function(lambda, lower, upper, where = "marginal") {
  if (!is.null(lower) && length(lower) != length(lambda)) {
    stop(sprintf("%s: lambda.lower must match length of coefficients (got %d vs %d).",
                 where, length(lower), length(lambda)), call. = FALSE)
  }
  if (!is.null(upper) && length(upper) != length(lambda)) {
    stop(sprintf("%s: lambda.upper must match length of coefficients (got %d vs %d).",
                 where, length(upper), length(lambda)), call. = FALSE)
  }
  if (!is.null(lower)) attr(lambda, "lower") <- lower
  if (!is.null(upper)) attr(lambda, "upper") <- upper
  lambda
}

