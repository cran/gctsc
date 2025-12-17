## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.width = 6,
  fig.height = 4
)
library(gctsc)
set.seed(1)

## -----------------------------------------------------------------------------
library(gctsc)

n  <- 100
mu <- 10
phi <- 0.2
arma_order <- c(1, 0)
tau <- c(phi)

# Simulate Poisson count data
sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 7)
y <- sim_data$y

plot(y, type = "l", main = "Simulated Poisson AR(1) Counts")


## -----------------------------------------------------------------------------
n <- 100
phi <- 0.5
tau <- c(phi)

beta_0 <- 1.2
prob   <- plogis(beta_0)
rho    <- 0.1
pi0    <- 0.8
size   <- 24

set.seed(7)
sim_data <- sim_zibb(prob = prob * rep(1,n),
                     rho = rho, pi0 = pi0,
                     size = size, tau = tau,
                     arma_order = c(1,0),
                     nsim = n)
y <- sim_data$y

plot(y, type = "l", main = "Simulated ZIBB AR(1) Counts")


## -----------------------------------------------------------------------------
n <- 300
mu <- 10
phi <- 0.5
theta <- 0.2
arma_order <- c(1,1)
tau <- c(phi, theta)

# Simulate Poisson count data
sim_data <- sim_poisson(mu, tau, arma_order, nsim = n, seed = 1)
y <- sim_data$y

# Compute bounds for truncated MVN using marginal object
marg <- poisson.marg()
X    <- matrix(1, nrow = n)
ab   <- marg$bounds(y, X, lambda = mu)
lower <- ab[,1]
upper <- ab[,2]

# Likelihood approximation
llk_tmet <- pmvn_tmet(lower, upper, tau, od = arma_order, pm = 30, QMC = TRUE)
llk_ghk  <- pmvn_ghk(lower, upper, tau, od = arma_order, QMC = TRUE)

c(TMET = llk_tmet, GHK = llk_ghk)


## -----------------------------------------------------------------------------
fit <- gctsc(
  formula  = y ~ 1,
  marginal = poisson.marg(),
  cormat   = arma.cormat(p = 1, q = 1),
  method   = "CE"
)

summary(fit)

predict(fit)


## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
pred_tmet <- predict(
  fit,
  method  = "TMET",
  y_obs   = 10
)

pred_tmet

## -----------------------------------------------------------------------------
## Load weekly Campylobacter incidence data
data("campyl", package = "gctsc")
y <- campyl[,"X1"]
n <- length(y)

## Plot the time series
ts_y <- ts(y, start = c(2001, 1), frequency = 52)
plot(ts_y, main = "Weekly Campylobacter Incidence",
     ylab = "Cases", xlab = "Year")

## ---------------------------------------------------------------
## Construct seasonal covariates
## ---------------------------------------------------------------
## Seasonal structure appears to have yearly periodicity,
## so we include sine/cosine terms with period = 52 weeks.

time <- 1:n
X <- data.frame(
  intercept = 1,
  sin52 = sin(2 * pi * time / 52),
  cos52 = cos(2 * pi * time / 52)
)

## Combine response and covariates
data_df <- data.frame(Y = y, X)

## Use the first 800 observations for model fitting
train_end <- 1000
data_train <- data_df[1:train_end, ]

## ---------------------------------------------------------------
## Fit a Negative Binomial Gaussian Copula model
## ---------------------------------------------------------------
## We use:
##   - Negative Binomial marginal (log link)
##   - ARMA(1,1) latent correlation structure
##   - GHK likelihood approximation
##
## The model is:
##     E[Y_t] = exp(β0 + β1 sin + β2 cos)
##
## Covariates enter only through the marginal mean.

fit <- gctsc(
  formula  = Y ~ sin52 + cos52,
  data     = data_train,
  marginal = negbin.marg(link = "log"),
  cormat   = arma.cormat(p = 1, q = 1),
  method   = "CE",  ### method can be changed to TMET or GHK
  options  = gctsc.opts(seed = 1)   # fixed seed for reproducibility
)

summary(fit)

## ---------------------------------------------------------------
## Residual diagnostics
## ---------------------------------------------------------------
plot(fit)

## ---------------------------------------------------------------
## One-step-ahead prediction
## ---------------------------------------------------------------
## Predict Y_{801} using fitted model
new_obs_index <- train_end + 1
pred <- predict(
  fit,
  X_test = data_df[new_obs_index, ],
  y_obs  = data_df[new_obs_index, "Y"]
)

pred


