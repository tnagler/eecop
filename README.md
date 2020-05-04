
# eecop

<!-- badges: start -->
<!-- badges: end -->

eecop is an R package for copula-based regression by the method of [Nagler and Vatter (2020)](https://arxiv.org/abs/1801.10576). A
model for marginal distributions and copula between response and covariates
is estimated. Predictions for quantiles or expectiles can then be derived
from solving a weighted estimating equations.

## Installation

You can install the released version of eecop from github with:

``` r
devtools::install_github("tnagler/eecop")
```

## Example

``` r
library(eecop)
## model with continuous variables
x <- matrix(rnorm(200), 100, 2)
y <- rowSums(x) + rnorm(100)

fit <- eecop(y, x)
predict(fit, x, t = c(0.5, 0.9), type = "quantile")
```

## References

Nagler, T. and Vatter, T. (2020). Solving estimating equations with copulas.
[*arXiv:1801.10576*](https://arxiv.org/abs/1801.10576)

