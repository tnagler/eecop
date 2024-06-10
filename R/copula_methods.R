#' @importFrom stats qnorm
#' @noRd
fit_w_normal <- function(v, u, weights) {
  rho <- wdm::wdm(qnorm(cbind(v, u)), weights = weights)
  c_YX <- copula::normalCopula(
    param = copula::P2p(rho),
    dim = ncol(cbind(v, u)),
    dispstr = "un"
  )
  c_Y_fix <- if (NCOL(v) > 1) {
    c_Y <- copula::normalCopula(
      param = copula::P2p(rho[1:NCOL(v), 1:NCOL(v)]),
      dim = NCOL(v),
      dispstr = "un"
    )
    copula::dCopula(v, c_Y)
  } else {
    rep(1, length(v))
  }

  function(u_new) {
    copula::dCopula(cbind(v, u_new), c_YX) / c_Y_fix
  }
}

#' @noRd
fit_w_vine <- function(v, u, weights, mult = 1, ...) {
  n <- NROW(u)
  d <- ncol(cbind(v, u))
  mult <- mult * (choose(d + 1, 2))^(1 / 6) * n^(-1 / 6 + 1 / 5)
  c_YX <- rvinecopulib::vinecop(cbind(v, u), weights = weights, mult = mult, ...)
  c_Y <- rvinecopulib::vinecop(v, weights = weights, mult = mult, ...)
  c_Y_fix <- rvinecopulib::dvinecop(v, c_Y)
  function(u_new) {
    rvinecopulib::dvinecop(cbind(v, u_new), c_YX) / c_Y_fix
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
  bws <- replicate(n, bw, simplify = FALSE)
  bws <- do.call(rbind, bws)

  if (length(weights) == 0) {
    weights <- rep(1, n)
  }

  if (d == 1) {
    f_Y_fix <- ks::dnorm.mixt(
      x = v,
      mus = v,
      sigmas = rep(sqrt(bws[1, 1]), n),
      props = weights / sum(weights)
    )
  } else {
    f_Y_fix <- ks::dmvnorm.mixt(
      x = v,
      mus = v,
      Sigmas = replicate(n, bws[1:d, 1:d], simplify = FALSE)
        |> do.call(rbind, args = _),
      props = weights / sum(weights)
    )
  }

  function(u_new) {
    f_YX <- ks::dmvnorm.mixt(
      x = cbind(v, u_new),
      mus = cbind(v, u),
      Sigmas = bws,
      props = weights / sum(weights)
    )
    f_YX / f_Y_fix
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
