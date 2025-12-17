## -------------------------------
## Example: Zero-Inflated Binomial AR(1) model
## -------------------------------
## Model parametrization note:
## ---------------------------
## For simulation, the user specifies:
##   - prob(t): Binomial success probability
##   - pi0(t):  zero–inflation probability
## These are probabilities on the (0,1) scale.

## For estimation, gctsc() uses a regression formulation:
##
##   logit( prob(t) ) = η_prob(t) = β_0 
##   logit( pi0(t)  ) = η_pi0(t)  = α_0 + α_1 z_1(t) + ...
##
## where:
##   prob(t) = plogis(η_prob(t))
##   pi0(t)  = plogis(η_pi0(t))
##
## This parametrization:
##   • ensures 0 < prob(t), pi0(t) < 1,
##   • allows covariates naturally,
##   • matches the GLM framework,
##   • is the scale used for MLE in gctsc().
##
## Thus, simulation uses probability parameters directly,
## while estimation uses their logistic regression
## representations through the model formula.



## -------------------------------
## Example: Zero-Inflated Binomial AR(1) 
## -------------------------------

library(gctsc)

## Model parametrization note:
## ---------------------------
## Simulation uses:
##   prob(t) = constant probability of success
##   pi0(t)  = constant zero-inflation probability
##
## Estimation uses logistic regression:
##   logit(prob(t)) = β_0          (intercept-only)
##   logit(pi0(t))  = α_0          (intercept-only)
##
## where prob(t) = plogis(β_0),  pi0(t) = plogis(α_0).

## --- Parameter setup ---
n    <- 500
size <- 24
prob <- 0.2           # Binomial success probability
pi0  <- 0.2           # Zero-inflation probability
phi  <- 0.8           # AR(1) parameter
tau  <- c(phi)
arma_order <- c(1, 0)

## --- Simulate ZIB count time series ---
set.seed(1)
sim_data <- sim_zib(
  prob       = rep(prob, n),
  pi0        = rep(pi0, n),
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  nsim       = n
)
y <- sim_data$y

## --- Fit Gaussian copula ZIB model using TMET ---
fit_zib <- gctsc(
  formula  = list(
    mu  = y ~ 1,   # logit(prob(t)) = β_0
    pi0 = ~ 1      # logit(pi0(t))  = α_0
  ),
  marginal = zib.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zib)
plot(fit_zib)
predict(fit_zib)


## -------------------------------
## Example: Zero-Inflated Binomial AR(1) with covariates
## -------------------------------

library(gctsc)

## Model parametrization note:
## ---------------------------
## Simulation uses a time-varying π0(t) directly on (0,1).
##
## Estimation uses logistic regression:
##   logit(pi0(t)) = α_0 + α_1 * Spring + α_2 * Summer + α_3 * Winter
##
## This allows π0(t) to vary seasonally through covariates.

n    <- 500
size <- 24
prob <- 0.2
phi  <- 0.8
tau  <- c(phi)
arma_order <- c(1, 0)

## --- Construct seasonal covariates for π0(t) ---
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))

X_pi <- model.matrix(~ season)
colnames(X_pi) <- make.names(colnames(X_pi))

## True zero-inflation function
beta_pi  <- c(0.2, 1, 0.3, 0.5)
logit_pi <- X_pi %*% beta_pi
pi0      <- plogis(logit_pi)

## --- Simulate ZIB time series ---
set.seed(1)
sim_data <- sim_zib(
  prob       = rep(prob, n),
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  nsim       = n
)
y <- sim_data$y

df <- data.frame(y = y, X_pi)

## --- Fit Gaussian copula ZIB model (TMET) ---
fit_zib_cov <- gctsc(
  formula  = list(
    mu  = y ~ 1,
    pi0 = ~ seasonSpring + seasonSummer + seasonWinter
  ),
  data     = df,
  marginal = zib.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zib_cov)
plot(fit_zib_cov)
predict(fit_zib_cov, X_test = df[200, ])
