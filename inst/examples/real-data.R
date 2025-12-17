## -------------------------------------------
## Example: Real Data – ZIBB Gaussian Copula Model
## Kickapoo Downtown Airport “Hot Hours”
## -------------------------------------------

library(gctsc)

## --- Load data ---
data("KCWC", package = "copTSC")
KCWC$date <- as.Date(KCWC$date)
y <- KCWC$hot

n <- length(y)
time <- 1:n

## Seasonal covariates for π0(t)
X_pi <- cbind(
  Intercept = 1,
  x_sin = sin(2 * pi * time / 365),
  x_cos = cos(2 * pi * time / 365)
)

## --- Train/Test split ---
n_train <- 500
y_train <- y[1:n_train]
X_pi_train <- X_pi[1:n_train, , drop = FALSE]

train_data <- data.frame(cbind(y_train, X_pi_train))

## ===========================================
##   Fit ZIBB Marginal + AR(1) Copula (TMET)
## ===========================================
fit_tmet <- gctsc(
  formula  = list(mu = y_train ~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zibb.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_tmet)
plot(fit_tmet)         # PIT, residuals, ACF plot

## --- One-step-ahead prediction ---
t_pred <- n_train + 1
pred_tmet <- predict(
  fit_tmet, 
  method = "TMET",
  y_obs  = y[t_pred],
  X_test = X_pi[t_pred, ]
)

pred_tmet

## ===========================================
##   Fit ZIBB Marginal (GHK for comparison)
## ===========================================
fit_ghk <- gctsc(
  formula  = list(mu = y_train ~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zibb.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  options  = gctsc.opts(seed = 1, M = 1000)
)

pred_ghk <- predict(
  fit_ghk, 
  method = "GHK",
  y_obs  = y[t_pred],
  X_test =  X_pi[t_pred, ]
)

pred_ghk

## ===========================================
##   Optional: Fit Zero-Inflated Binomial (ZIB)
## ===========================================
fit_zib <- gctsc(
  formula  = list(mu = y_train ~ 1, pi0 = ~ x_sin + x_cos),
  data     = train_data,
  marginal = zib.marg(link = "logit", size = 24),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zib)
