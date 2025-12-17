#' @importFrom stats model.frame model.response model.matrix terms na.pass
#' @importFrom stats runif optim
#' @importFrom stats glm glm.fit binomial poisson fitted coef
#' @importFrom stats dbinom pbinom dpois ppois dnbinom pnbinom
#' @importFrom stats plogis qlogis as.formula
#' @importFrom stats acf pacf residuals sd toeplitz
#' @importFrom stats qnorm pnorm rnorm qbinom qpois qnbinom
#' @importFrom graphics par abline mtext hist
#' @importFrom grDevices dev.interactive
#' @importFrom utils data
#' @importFrom stats arima printCoefmat
#' @importFrom stats quasibinomial setNames
#' @importFrom utils tail
utils::globalVariables(c("lw", "n"))
NULL
