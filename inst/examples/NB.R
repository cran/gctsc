## -------------------------------
## Example: negative binomial ARMA(1,1) model
## -------------------------------

## --- Parameter setup ---
n <- 500
mu <- 1
phi <- 0.5
theta <- 0.8
arma_order <- c(1, 1)
tau <- c(phi,theta)
dispersion <- 2

## --- Simulate data ---
set.seed(7)
sim_data <- sim_negbin(mu = mu * rep(1,n), dispersion,
                        tau = tau,
                        arma_order = arma_order,
                        nsim = n)
y <- sim_data$y
X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- negbin.marg()
ab <- marg$bounds(y, X, c(mu,dispersion))

## --- Likelihood approximation ---
llk_tmet <- pmvn_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      pm = 30, QMC = TRUE)

llk_ghk  <- pmvn_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE)

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit Gaussian copula model using TMET ---
fit_tmet <- gctsc(
  formula = y ~ 1,
  marginal = negbin.marg(lambda.lower = c(0,0)),
  cormat   = arma.cormat(p = 1, q = 1),
  method   = "TMET",
  QMC      = TRUE,
  options = gctsc.opts(M = 1000, seed = 42)
)
summary(fit_tmet)
plot(fit_tmet)       # residual diagnostics
predict(fit_tmet)    # one-step forecasting

## --- Fit Gaussian copula model using GHK ---
fit_ghk <- gctsc(
  formula = y ~ 1,
  marginal = negbin.marg(lambda.lower = c(0,0)),
  cormat   = arma.cormat(p = 1, q = 1),
  method   = "GHK",
  QMC      = TRUE,
  options = gctsc.opts(M = 1000, seed = 42)
)

plot(fit_ghk)



## -------------------------------
## Example: Negative Binomial AR(1) model with covariates
## -------------------------------

library(gctsc)

n <- 300
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

## Covariate regression coefficients
beta <- c(1, 0.3, 1, 0.5)

## Negative binomial dispersion (gctsc parameterization: variance = mu + mu^2/dispersion)
dispersion <- 2

## Generate covariates (seasonal + autoregressive component)
set.seed(1)
zeta <- rnorm(n)
xi <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
}

X <- as.matrix(data.frame(
  x1 = rep(1, n),
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))

## Compute mean function
mu <- as.vector(exp(X %*% beta))

## Simulate Negative Binomial counts with AR(1) latent process
sim_data <- sim_negbin(
  mu        = mu,
  dispersion = dispersion,
  tau       = tau,
  arma_order = arma_order,
  nsim      = n,
  seed      = 10
)
y <- sim_data$y

## Assemble data for model fitting
data_df <- data.frame(Y = y, X)

## Fit Gaussian Copula model using GHK
fit_nb <- gctsc(
  formula  = Y ~ x2 + x3 + x4,
  data     = data_df[1:499,],
  marginal = negbin.marg(link = "log"),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_nb)
plot(fit_nb)

## One-step-ahead prediction
predict(fit_nb,y =data_df$Y[500],  X_test = data_df[500, c("x2","x3","x4")])
