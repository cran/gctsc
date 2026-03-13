## -------------------------------
## Example: Beta–Binomial AR(1) model
## -------------------------------

library(gctsc)

## --- Parametrization note ----------------------------------------------
## Simulation:
##   sim_bbinom() uses:
##       prob(t) ∈ (0,1)          (success probability)
##       rho ∈ (0,1)              (intra-class correlation)
##
## Estimation:
##   gctsc() fits the Beta–Binomial marginal using a logit link:
##       logit{prob(t)} = η_prob(t)
##   so that prob(t) = plogis(η_prob(t)).
##
## This parametrization:
##   - matches GLM-style modeling for Beta–Binomial mean structure,
##   - allows covariates to enter prob(t) naturally,
##   - prevents boundary issues with prob ∈ (0,1),
##   - keeps rho fixed or treated as a separate dispersion parameter.

## --- Parameter setup ---
n    <- 500
size <- 24
beta0 <- 0.2                  # logit-scale intercept
prob <- plogis(beta0)         # simulation-scale probability
rho  <- 0.18                  # intra-class correlation
phi  <- 0.8                   # AR(1) parameter
tau  <- c(phi)
arma_order <- c(1, 0)
df <- 5

## --- Simulate Beta–Binomial time series ---
set.seed(1)
sim_data <- sim_bbinom(
  prob       = rep(prob, n),
  rho        = rho,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "t",
  df =5,
  nsim       = n
)
y <- sim_data$y


X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- bbinom.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, c(prob,rho),family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit Gaussian copula Beta–Binomial model using TMET ---
fit_bbinom <- gctsc(
  formula  = y ~ 1,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET", family = "t", df=5,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_bbinom)
plot(fit_bbinom)
predict(fit_bbinom)

## -------------------------------
## Example: Beta–Binomial AR(1) with covariates
## -------------------------------

library(gctsc)

## --- Parametrization note ----------------------------------------------
## prob(t) is generated from a logistic regression:
##       logit{prob(t)} = X(t)^T β
## prob(t) = plogis(Xβ)
##
## The same logit parametrization is used in gctsc():
##   marginal = bbinom.marg(link = "logit")
##
## This ensures:
##   - directly comparable coefficients (β vs fitted β̂),
##   - natural covariate inclusion,
##   - probability always within (0,1).

## --- Parameter setup ---
n    <- 500
size <- 24
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)

## Overdispersion / intra-class correlation parameterization
## rho_vgam corresponds to 1 / (1 + θ) in VGAM
rho_vgam <- 1 / (1 + 5)

## --- Construct covariates (seasonal + AR structure) ---
zeta <- rnorm(n)
xi   <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j-1] - 0.4 * xi[j-2] + zeta[j]
}

X <- data.frame(
  x1 = 1,
  x2 = sin(2*pi*(1:n)/12),
  x3 = cos(2*pi*(1:n)/12),
  x4 = xi
)

## True logit(prob)
beta_true <- c(0.2, 0.3, 0.5, 0.3)
logit_prob <- as.matrix(X) %*% beta_true
prob <- plogis(logit_prob)

## --- Simulate Beta–Binomial time series ---
set.seed(1)
sim_data <- sim_bbinom(
  prob       = prob,
  rho        = rho_vgam,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "gaussian",
  nsim       = n
)
y <- sim_data$y



## --- Compute truncation bounds ---
marg <- bbinom.marg(link = "logit", size= size)
ab <- marg$bounds(y, as.matrix(X), c(0.2, 0.3, 0.5, 0.3, rho),family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)


## --- Fit t copula Beta–Binomial model ---
data_df <- data.frame(y = y, X)

fit_bbinom_cov <- gctsc(
  formula  = y ~ x2 + x3 + x4,
  data     = data_df,
  marginal = bbinom.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "t", df=5,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_bbinom_cov)
plot(fit_bbinom_cov)
predict(fit_bbinom_cov, X_test = data_df[200, ])

