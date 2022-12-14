
# BVshear

<!-- badges: start -->
<!-- badges: end -->

The goal of BVshear is to provide density, probability, quantile and random generation functions for the Brutsaert-Veron shear driven distribution.

## Installation

You can install the development version of BVshear from [GitHub](https://github.com/eliardocosta) with:

``` r
install.packages("remotes")
library(remotes)
install_github("eliardocosta/BVshear")
```

## Example

These are basic examples showings how to use the functions of the package for the Brutsaert-Veron shear driven distribution:

``` r
library(BVshear)
# density at 0.03 for u = 0.5
dBVshear(x = 0.03, u = 0.5) 

# log-density at 0.03 for u = 0.5
dBVshear(x = 0.03, u = 0.5, log = TRUE) 

# P(X <= 0.03) for u = 0.5
pBVshear(q = 0.03, u = 0.5) 

# P(X > 0.03) for u = 0.5
pBVshear(q = 0.03, u = 0.5, lower.tail = FALSE) 

# median of the distribution for u = 0.5
qBVshear(p = 0.5, u = 0.5)

set.seed(1234)
tau <- rBVshear(n = 1E4, u = 0.5)
eps <- 0.005
hist(tau, freq = FALSE, breaks = seq(0, max(tau) + eps, eps), xlim = c(0, 0.2), 
     ylim = c(0, 40), xlab = "Contact time", main = "Histogram of contact times")
```

