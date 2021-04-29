#' @importFrom stats qnorm
#' @noRd
fit_w_normal <- function(v, u, weights) {
  rho <- wdm::wdm(qnorm(cbind(v, u)), weights = weights)
  c_YX <- copula::normalCopula(
    param = copula::P2p(rho),
    dim = ncol(cbind(v, u)),
    dispstr = "un"
  )
  c_Y <- if (NCOL(v) > 1) {
    copula::normalCopula(
      param = copula::P2p(rho[1:NCOL(v), 1:NCOL(v)]),
      dim = ncol(v),
      dispstr = "un"
    )
  } else {
    function(v) rep(1, length(v))
  }

  function(u_new) {
    copula::dCopula(cbind(v, u_new), c_YX) / copula::dCopula(u_new, c_Y)
  }
}

#' @noRd
fit_w_vine <- function(v, u, weights, mult = 1, ...) {
  n <- NROW(u)
  mult <- mult * n^(1 / 5 - 1 / 4.5)
  c_YX <- rvinecopulib::vinecop(cbind(v, u), weights = weights, mult = mult, ...)
  c_Y <- rvinecopulib::vinecop(v, weights = weights, mult = mult, ...)
  function(u_new) {
    rvinecopulib::dvinecop(cbind(v, u_new), c_YX) / rvinecopulib::dvinecop(v, c_Y)
  }
}

#' @importFrom stats cov
fit_w_kde <- function(v, u, weights, mult = 1, ...) {
  u <- as.matrix(u)
  v <- as.matrix(v)
  n <- nrow(u)
  d <- ncol(v)
  p <- ncol(u)
  bw <- mult * n^(-2 / (4 + p)) * cov(cbind(v, u))
  function(u_new) {
    if (length(weights) == 0) {
      weights <- rep(1, n)
    }
    bws <- replicate(n, bw, simplify = FALSE)
    bws <- do.call(rbind, bws)
    f_YX <- ks::dmvnorm.mixt(
      x = u_new,
      mus = u,
      Sigmas = bws,
      props = weights / sum(weights)
    )
    if (d == 1) {
      f_Y <- ks::dnorm.mixt(
        x = v,
        mus = v,
        sigmas = rep(sqrt(bws), n),
        props = weights / sum(weights)
      )
    } else {
      f_Y <- ks::dmvnorm.mixt(
        x = v,
        mus = v,
        Sigmas = bws[1:d, 1:d],
        props = weights / sum(weights)
      )
    }

    f_YX / f_Y
  }
}


# fit_w_x_trafo_kde <- function(v, u, weights, mult = 1, q = 1, ...) {
# ## incomplete
#   y <- qnorm
#   x <- qnorm(u)
#   n <- nrow(x)
#   d <- ncol(x)
#   bw <- mult * n^(-2 / (4 + d)) * cov(x)
#   function(u_new) {
#     x_new <- qnorm(u_new)
#     n <- nrow(x)
#     if (length(weights) == 0) {
#       weights <- rep(1, n)
#     }
#     bws <- replicate(n, bw, simplify = FALSE)
#     bws <- do.call(rbind, bws)
#     f <- ks::dmvnorm.mixt(
#       x = x_new,
#       mus = x,
#       Sigmas = bws,
#       props = weights / sum(weights)
#     )
#     f / pmax(1e-20, apply(dnorm(x_new), 1, prod))
#   }
# }
