## -------------------------------
## Example: Zero-Inflated Beta-Binomial AR(1) model
## -------------------------------

library(gctsc)

## --- Parametrization note ----------------------------------------------
## Simulation:
##   sim_zibb() uses:
##       prob(t) ∈ (0,1)        (success probability for BB component)
##       rho ∈ (0,1)            (intra-class correlation)
##       pi0(t) ∈ (0,1)         (zero inflation)
##
## Estimation:
##   gctsc() fits two logistic regressions:
##       logit{prob(t)} = η_mu(t)
##       logit{pi0(t)}  = η_pi0(t)
##
## allowing covariates in BOTH the mean and zero-inflation components.
## This avoids boundary problems and ensures prob(t), pi0(t) ∈ (0,1).

## --- Parameter setup ---
n    <- 500
size <- 24
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)

## True parameters
beta_mu  <- 1.2                 # logit-scale parameter for prob
prob     <- plogis(beta_mu)     # simulation-scale success probability
rho      <- 0.1                 # BB overdispersion (ICC)
pi0      <- 0.2                 # constant zero-inflation prob
df= 10

## --- Simulate ZIBB time series ---
set.seed(7)
sim_data <- sim_zibb(
  prob       = rep(prob, n),
  rho        = rho,
  pi0        = rep(pi0, n),
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "t",
  df= 10,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = matrix(1, nrow = n), pi0 = matrix(1, nrow = n))

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), qlogis(pi0))

marg <- zib.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)
## --- Fit ZIBB copula model using TMET ---
fit_zibb <- gctsc(
  formula  = list(mu = y ~ 1, pi0 = ~ 1),
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK", family = "t", df= 10,
  options  = gctsc.opts(seed = 7, M = 1000)
)

summary(fit_zibb)
plot(fit_zibb)
predict(fit_zibb)

## -------------------------------
## Example: ZIBB AR(1) with seasonal zero-inflation π₀(t)
## -------------------------------
## --- Parametrization note (for Example 2) ------------------------------------
## Simulation:
##   prob(t)  is constant here (no covariates in μ), passed directly to sim_zibb().
##   pi0(t)   varies with seasonal covariates via:
##        logit{pi0(t)} = X_pi(t)^T β_pi
##   rho      controls Beta-Binomial overdispersion.

## Estimation:
##   gctsc() fits:
##        logit{prob(t)} = η_mu(t)
##        logit{pi0(t)}  = η_pi0(t)
##   even if prob(t) is constant in the simulation.
##
## This ensures:
##   – valid probabilities prob(t), pi0(t) ∈ (0,1)
##   – natural support for covariates in the zero-inflation model
##   – consistency across all ZI marginals in the package.
library(gctsc)

## --- Parameter setup ---
n    <- 500
size <- 24
phi  <- 0.85
tau  <- c(phi)
arma_order <- c(1, 0)
df= 10

prob <- plogis(0.2)        # constant prob(t)
rho  <- 0.16               # ICC for BB component

## Seasonal covariates for π₀(t)
X_pi <- cbind(
  1,
  sin(2*pi*(1:n)/12),
  cos(2*pi*(1:n)/12)
)
colnames(X_pi) <- c("int", "sin", "cos")

beta_pi <- c(1, 2, 5)
pi0 <- plogis(X_pi %*% beta_pi)

## --- Simulate ZIBB data ---
set.seed(1)
sim_data <- sim_zibb(
  prob       = rep(prob, n),
  rho        = rho,
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  family = "t",
  df= 10,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = matrix(1, nrow = n), pi0 = as.matrix(X_pi))

## --- Compute truncation bounds ---
lambda <- c(qlogis(prob), qlogis(rho), beta_pi)

marg <- zibb.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)
## --- Fit using formula interface ---
df <- data.frame(y = y, X_pi)

fit_zibb_seasonal <- gctsc(
  formula = list(mu = y ~ 1,
                 pi0 = ~ sin + cos),
  data     = df,
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "t", df= 10,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zibb_seasonal)
plot(fit_zibb_seasonal)
predict(fit_zibb_seasonal, X_test = df[200, ])


## -------------------------------
## Example: ZIBB AR(1) with covariates in μ(t) and π₀(t)
## -------------------------------
## --- Parametrization note (for Example 3) ------------------------------------
## Simulation:
##   sim_zibb() uses:
##      prob(t) = plogis( X_mu(t)^T β_mu )        # mean function on probability scale
##      pi0(t)  = plogis( X_pi(t)^T α_pi )        # zero-inflation on probability scale
##      rho     is the Beta-Binomial ICC.

## Estimation:
##   gctsc() fits two logistic regressions:
##       logit{prob(t)} = X_mu(t)^T β_mu
##       logit{pi0(t)}  = X_pi(t)^T α_pi
##
## Notes:
##   – This automatically keeps prob(t), pi0(t) within (0,1).
##   – Allows arbitrary covariates in BOTH μ(t) and π₀(t).
##   – This parametrization is standard for zero-inflated count models.

library(gctsc)

## --- Parameter setup ---
n    <- 1500
size <- 10
phi  <- 0.5
tau  <- c(phi)
arma_order <- c(1, 0)
rho  <- 0.3
df= 10

## --- Covariates for μ(t) ---
X_mu <- cbind(
  1,
  sin(2*pi*(1:n)/365),
  cos(2*pi*(1:n)/365)
)
colnames(X_mu) <- c("int", "sin", "cos")

beta_mu <- c(1.2, 2, 3)
prob <- plogis(X_mu %*% beta_mu)

## --- Covariates for π₀(t) ---
day_of_year <- rep(1:365, length.out = n)
season <- factor(ifelse(day_of_year < 100, "Winter",
                        ifelse(day_of_year < 180, "Spring",
                               ifelse(day_of_year < 270, "Summer", "Fall"))))
X_pi <- model.matrix(~ season)

alpha_pi <- c(0.2, 0.3, 0.5, 0.7)
pi0 <- plogis(X_pi %*% alpha_pi)

## --- Simulate ZIBB data ---
set.seed(1)
sim_data <- sim_zibb(
  prob       = prob,
  rho        = rho,
  pi0        = pi0,
  size       = size,
  tau        = tau,
  arma_order = arma_order,
  nsim       = n
)
y <- sim_data$y

X <- list(mu = as.matrix(X_mu), pi0 = as.matrix(X_pi))

## --- Compute truncation bounds ---
lambda <- c(beta_mu, qlogis(rho), alpha_pi)

marg <- zibb.marg(link = "logit", size= size)
ab <- marg$bounds(y, X, lambda, family ="t", df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order, 
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df= df)

c(TMET = llk_tmet, GHK = llk_ghk)


## --- Fit ZIBB copula model using formula interface ---
data_df <- data.frame(y = y, X_mu, X_pi)

fit_zibb_cov <- gctsc(
  formula = list(mu  = y ~ sin + cos,
                 pi0 = ~ seasonSpring + seasonSummer + seasonWinter),
  data     = data_df,
  marginal = zibb.marg(link = "logit", size = size),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "TMET",family = "t", df= 10,
  options  = gctsc.opts(seed = 1, M = 1000)
)

summary(fit_zibb_cov)
plot(fit_zibb_cov)

predict(fit_zibb_cov, X_test = data_df[500, ])

