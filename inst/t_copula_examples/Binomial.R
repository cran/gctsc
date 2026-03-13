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
df <- 10

X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- binom.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, prob,family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Simulate Binomial count time series ---
set.seed(1)
sim_data <- sim_binom(
  prob       = rep(prob, n),     # simulation scale
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family     = "t",
  df         = df,
  nsim       = n
)
y <- sim_data$y

## --- Fit Student t copula Binomial model using GHK ---
fit_binom <- gctsc(
  formula  = y ~ 1,
  marginal = binom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),family = "t", df= 10,
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_binom)
plot(fit_binom)
predict(fit_binom)
