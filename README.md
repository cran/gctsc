---
output:
  pdf_document: default
  html_document: default
---
# gctsc

`gctsc` provides fast and scalable likelihood inference for **Gaussian and Student–t copula models for count time series**.

The package supports a wide range of discrete marginals:

- Poisson  
- Negative Binomial  
- Binomial  
- Beta Binomial  
- Zero-Inflated Poisson (ZIP)  
- Zero-Inflated Binomial (ZIB)  
- Zero-Inflated Beta-Binomial (ZIBB)

The latent dependence structure is modeled via **ARMA(p, q)** processes.

Likelihood evaluation is performed using one of the following approximation methods:

- **TMET** — Time Series Minimax Exponential Tilting  
- **GHK** — Geweke–Hajivassiliou–Keane simulator  
- **CE** — Continuous Extension  

The implementation exploits ARMA structure for efficient high-dimensional computation.

Additional features include:

- Simulation utilities  
- Randomized quantile residual diagnostics  
- Predictive distributions and scoring rules (CRPS, LOGS)  

---

## Installation

From CRAN (after release):

```r
install.packages("gctsc"
```

From Github:
remotes::install_github("QNNHU/gctsc")

## Quick Example:
```r
library(gctsc)

# Simulate Poisson AR(1) data under a Gaussian copula
set.seed(1)
y <- sim_poisson(
  mu = 5,
  tau = 0.5,
  arma_order = c(1, 0),
  nsim = 300,
  family = "gaussian"
)$y

# Fit model
fit <- gctsc(
  y ~ 1,
  data = data.frame(y = y),
  marginal = poisson.marg(),
  cormat = arma.cormat(p = 1, q = 0),
  method = "TMET",
  family = "gaussian",
  options = gctsc.opts(M = 1000)
)

summary(fit)

# Diagnostic plots
plot(fit)

# One-step prediction
predict(fit)
```
# What Makes gctsc Different?

Compared to existing implementations, `gctsc` added:

- Exploits ARMA structure for scalable likelihood evaluation in time series settings

- Supports zero-inflated marginals with flexible covariate specification, including seasonal components

- Implements scalable minimax exponential tilting (TMET) for efficient likelihood approximation

- Provides a linear-cost GHK importance sampling implementation

- Implements fast continuous extension method

- Supports Student–t copulas for modeling heavy-tailed dependence

- Computes full predictive distributions for discrete time series

# References

If you use this package in published work, please cite:

Nguyen, Q. N., & De Oliveira, V. (2026).
Approximating Gaussian copula models for count time series: Connecting the distributional transform and a continuous extension.
Journal of Applied Statistics.

Nguyen, Q. N., & De Oliveira, V. (2026).
Likelihood Inference in Gaussian Copula Models for Count Time Series via Minimax Exponential Tilting.
Computational Statistics & Data Analysis.

Nguyen, Q. N., & De Oliveira, V. (2026).
Scalable Likelihood Inference for Student–t Copula Count Time Series.
Manuscript in preparation.
