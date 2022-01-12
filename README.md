
<!-- README.md is generated from README.Rmd. Please edit that file -->

# flatlandr

<!-- badges: start -->
<!-- badges: end -->

This package provides functions for chebyshev interpolation in 2
dimensions

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mdgordo/flatlandr")
```

## Example

``` r
library(flatlandr)

## function for interpolation:
fxy <- function(x, y) {
  r = log(x) + x^2*y - y^3
}

### create 10 nodes in each dimension on (0,10) and (5,15)
lb = c(0,5)
ub = c(10,15)
ss = chebstatespace(10, lb, ub)

### apply function
y = mapply(fxy, ss$x, ss$y)
y = matrix(y, ncol = 10)

### find coefficient matrix
cmat = chebcoefs(y, degree = 5)

### interpolate
chebpred(2.5, 10.3, cmat, lb, ub)
```