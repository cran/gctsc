## -------------------------------
## Example: Poisson AR(1) model
## -------------------------------

## --- Parameter setup ---
n <- 1000
mu <- 10
phi <- 0.2
arma_order <- c(1, 0)
tau <- c(phi)
df <- 15
family <- "t"
## --- Simulate data ---
set.seed(7)
sim_data <- sim_poisson(mu = mu,
                        tau = tau,
                        arma_order = arma_order,
                        family = "t",
                        df = df,
                        nsim = n)
y <- sim_data$y
X <- matrix(1, nrow = n)

## --- Compute truncation bounds ---
marg <- poisson.marg()
ab <- marg$bounds(y, X, mu, family = family, df= df)

## --- Likelihood approximation ---
llk_tmet <- pmvt_tmet(lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      pm = 30, QMC = TRUE, df= df)

llk_ghk  <- pmvt_ghk( lower = ab[,1], upper = ab[,2],
                      tau = tau, od = arma_order,
                      QMC = TRUE, df=df)

c(TMET = llk_tmet, GHK = llk_ghk)

## --- Fit t copula model using TMET ---
fit_CE <- gctsc(
  formula = y ~ 1,
  marginal = poisson.marg(lambda.lower = 0),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "CE",
  family   = "t",
  df= 15,
  QMC      = TRUE
)

plot(fit_CE)       # residual diagnostics
predict(fit_CE)    # one-step forecasting

## --- Fit t copula model using GHK ---
fit_ghk <- gctsc(
  formula = y ~ 1,
  marginal = poisson.marg(lambda.lower = 0),
  cormat   = arma.cormat(p = 1, q = 0),
  method   = "GHK",
  family   = "t",
  df= 15,
  QMC      = TRUE
)

plot(fit_ghk)

## -------------------------------
## Example: Poisson AR(1) model with covariates
## -------------------------------

n <- 500
phi <- 0.8
tau <- c(phi)
arma_order <- c(1, 0)

## --- Generate covariates (seasonal + autoregressive) ---
set.seed(1)
zeta <- rnorm(n)
xi <- numeric(n)
for (j in 3:n) {
  xi[j] <- 0.6 * xi[j - 1] - 0.4 * xi[j - 2] + zeta[j]
}

X <- as.matrix(data.frame(
  x1 = 1,
  x2 = sin(2 * pi * (1:n) / 12),
  x3 = cos(2 * pi * (1:n) / 12),
  x4 = xi
))

beta <- c(0.1, 0.3, 1, 3)
mu <- exp(X %*% beta)

## --- Simulate Poisson response ---
sim_data <- sim_poisson(mu = mu, tau = tau,
                        arma_order = arma_order, family   = "t", df= df,
                        nsim = n, seed = 1)
y <- sim_data$y

## --- Compute bounds and log-likelihood approximations ---
marginal <- poisson.marg(link = "log")
ab <- marginal$bounds(y, X, beta, family   = "t", df= df)

llk_tmet_qmc <- pmvn_tmet(
  lower = ab[, 1],
  upper = ab[, 2],
  tau   = tau,
  od    = arma_order
)

## --- Fit t copula model ---
data_df <- data.frame(Y = y, X)
data_train <- data_df[1:400,]
fit <- gctsc(
  formula  = Y ~ x2 + x3 + x4,
  data     = data_train,
  marginal = poisson.marg(link = "log"),
  cormat   = arma.cormat(p = 1, q = 0),
  family   = "t",
  method   = "GHK",
  df= 15,
  options  = gctsc.opts(seed = 1)
)

summary(fit)
plot(fit)
predict(fit, X_test = ( data_df[401,3:5]), y_obs =  data_df[401,"Y"])
