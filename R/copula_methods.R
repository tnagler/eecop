#' @importFrom stats qnorm
#' @noRd
fit_copula_normal <- function(u, weights) {
  rho <- wdm::wdm(qnorm(u), weights = weights)
  cop_model <- copula::normalCopula(
    param = copula::P2p(rho),
    dim = ncol(u),
    dispstr = "un"
  )
  function(u_new) copula::dCopula(u_new, cop_model)
}

#' @noRd
fit_copula_vine <- function(u, weights, mult = 1, ...) {
  n <- NROW(u)
  mult <- mult * n^(1 / 5 - 1 / 4.5)
  cop_model <- rvinecopulib::vinecop(u, weights = weights, mult = mult, ...)
  function(u_new) rvinecopulib::dvinecop(u_new, cop_model)
}

#' @importFrom stats qnorm dnorm cov
#' @noRd
fit_copula_kde <- function(u, weights, mult = 1, q = 1, ...) {
  x <- qnorm(u)
  n <- nrow(x)
  d <- ncol(x)
  bw <- mult * n^(-2 / (4 + d)) * cov(x)
  function(u_new) {
    x_new <- qnorm(u_new)
    n <- nrow(x)
    if (length(weights) == 0) {
      weights <- rep(1, n)
    }
    bws <- replicate(n, bw, simplify = FALSE)
    bws <- do.call(rbind, bws)
    f <- ks::dmvnorm.mixt(
      x = x_new,
      mus = x,
      Sigmas = bws,
      props = weights / sum(weights)
    )
    f / max(1e-20, apply(dnorm(x_new), 1, prod))
  }
}
