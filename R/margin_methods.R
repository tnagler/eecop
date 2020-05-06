#' @importFrom stats qnorm pnorm sd
fit_margin_normal <- function(x, weights) {
  if (length(weights) == 0) {
    weights <- rep(1, length(x))
    m <- mean(x)
    s <- sd(x)
  } else {
    weights <- weights / sum(weights)
    m <- sum(x * weights)
    s <- sqrt(weights %*% (x - m)^2)
  }
  function(y) pnorm(y, m, s)
}

#' @importFrom kde1d kde1d pkde1d
fit_margin_kde <- function(x, weights) {
  mult <- length(x)^(1 / 5 - 1 / 3.5)
  fit <- kde1d::kde1d(x, weights = weights, mult = mult, deg = 0)
  function(y) kde1d::pkde1d(y, fit)
}
