
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gammaFrailty

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/gammaFrailty)](https://CRAN.R-project.org/package=gammaFrailty)
<!-- badges: end -->

The gammaFrailty package allows the estimation of serially correlated
Gamma-frailty models for count data as described in Henderson &
Shimakura (2003).

## Installation

You can install the development version of gammaFrailty via:

``` r
# install.packages("devtools")
devtools::install_github("giuseppealfonzetti/gammaFrailty")
```

## Example

Simple illustrative setting taken from Henderson & Shimakura (2003)

``` r
library(gammaFrailty)
seed <- 123
set.seed(seed)

# Setup data generating process
## Set dimensions
p <- 12; q <- 8; m <- 2; n <- 250

## Set parameters
xi <- 2/q; rho <- .6;
int <- rep(log(4), p)
b <- rep(0, m) 

## True parameter vector (and its reparameterisation)
theta <- c(xi, rho, b, int)
repar_theta <- partorepar(theta)

## Generate binary covariates
X <- matrix(rbinom(m*n, 1, .5), n, m)

## Generate data
data <- generate_data(
    INTERCEPT = int,
    BETA = b,
    X = X,
    Q = q,
    RHO = rho,
    SEED = seed
)


# Estimation
## Initialisation vector
par_init <- repar_theta + runif(length(repar_theta), -1, 1)

fit <- fit_gammaFrailty(
    DATA_LIST = list('DATA' = data, 'X' = X),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    VERBOSEFLAG= 0,
    INIT = par_init,
    ITERATIONS_SUBSET = NULL
)
#> 1. Initialising at init vector.
#> 2. Optimising with ucminf...
#> 3. Done! (0.94 secs)

mean((repartopar(fit$theta)-theta)^2)
#> [1] 0.001821759
mean((repartopar(par_init)-theta)^2)
#> [1] 0.3043935
```
