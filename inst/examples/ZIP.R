## -------------------------------
## Example: Zero-Inflated Poisson AR(1) model
## -------------------------------

library(gctsc)

## --- Parametrization note ----------------------------------------------
## Simulation:
##    mu        = Poisson mean (identity scale)
##    pi0(t)    = zero-inflation probability in (0,1)
##
## Estimation (in gctsc):
##    logit{pi0(t)} = η_pi0(t) = α_0 + α_1 x_1(t) + ...
## so that:
##    pi0(t) = plogis(η_pi0(t))
##
## This allows covariates and ensures pi0(t) ∈ (0,1).

## --- Parameter setup ---
n   <- 500
mu  <- 10               # Poisson mean
pi0 <- 0.4              # zero-inflation probability (probability scale)
phi <- 0.5
tau <- c(phi)
arma_order <- c(1, 0)

## --- Simulate ZIP data (pi0 on probability scale) ---
set.seed(1)
sim_data <- sim_zip(
  mu         = rep(mu, n),
  pi0        = rep(pi0, n),
  tau        = tau,
  arma_order = arma_order,
  nsim       = n
)
y <- sim_data$y

## --- Fit ZIP Gaussian copula model using TMET ---
## zip.marg(link="identity") means:
##   μ(t)   = identity(X_mu β)
##   logit(pi0(t)) = η_pi0(t)
## Estimation always uses logit for pi0.
fit_zip <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ 1),
  marginal = zip.marg(link = "identity"),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "CE",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zip)
plot(fit_zip)
predict(fit_zip)
