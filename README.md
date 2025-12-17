gctsc
================

# gctsc

`gctsc` provides fast and scalable likelihood inference for **Gaussian
copula models for count time series**, supporting a wide range of
marginals:

- Poisson  
- Negative Binomial  
- Binomial  
- Beta Binomial
- Zero-Inflated Poisson (ZIP)  
- Zero-Inflated Binomial (ZIB)  
- Zero-Inflated Beta-Binomial (ZIBB)

The package implements several likelihood approximation methods —
including the proposed  
**TMET** (Time Series Minimax Exponential Tilting) and **GHK** — and
exploits the **ARMA dependence structure** for efficient
high-dimensional computation.

Additional features include simulation utilities, residual diagnostic
tools, and one-step prediction.

------------------------------------------------------------------------

## Reference

If you use this package in published work, please cite:

> Nguyen, N. & De Oliveira, V. (2025).  
> *Likelihood Inference in Gaussian Copula Models for Count Time Series
> via Minimax Exponential Tilting.*
