## -------------------------------
## Example: Binomial AR(1) model
## -------------------------------

library(gctsc)

## --- Parametrization note ----------------------------------------------
## Simulation:
##   prob(t) ∈ (0,1) is used directly in sim_binom()
##
## Estimation:
##   gctsc() fits the Binomial marginal using a logit link:
##       logit{prob(t)} = η_prob(t)
##   so that:
##       prob(t) = plogis(η_prob(t))
##
## This parametrization:
##   - allows covariates to enter prob(t) naturally,
##   - avoids boundary issues (prob cannot hit 0 or 1),
##   - matches standard GLM/binomial practice.

## --- Parameter setup ---
n    <- 200
size <- 24                # number of trials
prob <- 0.3               # success probability (probability scale)
phi  <- 0.8               # AR(1) dependence
tau  <- c(phi)
arma_order <- c(1, 0)

## --- Simulate Binomial count time series ---
set.seed(1)
sim_data <- sim_binom(
  prob       = rep(prob, n),     # simulation scale
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "gaussian",
  nsim       = n
)
y <- sim_data$y

## --- Fit Gaussian copula Binomial model using GHK ---
fit_binom <- gctsc(
  formula  = y ~ 1,
  marginal = binom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),family = "gaussian",
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_binom)
plot(fit_binom)
predict(fit_binom)
