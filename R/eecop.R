#' Title
#'
#' @param y
#' @param x
#' @param copula_method
#' @param margin_method
#' @param weights
#' @param ...
#'
#' @export
#'
#' @importFrom assertthat assert_that is.string
#'
#' @examples
eecop <- function(y, x, copula_method = "vine", margin_method = "kde",
                  weights = numeric(), ...) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  assert_that(
    nrow(y) == nrow(x),
    is.string(copula_method),
    is.string(margin_method),
    copula_method %in% c("vine", "normal", "kde", "bernstein"),
    margin_method %in% c("kde", "normal"),
    is.numeric(weights)
  )

  n <- nrow(y)
  q <- ncol(y)
  p <- ncol(x)

  margins_Y <- lapply(
    seq_len(q),
    function(j) fit_margin(y[, j], margin_method, weights)
  )
  margins_X <- lapply(
    seq_len(p),
    function(j) fit_margin(x[, j], margin_method, weights)
  )

  V <- sapply(seq_len(q), function(j) margins_Y[[j]](y[, j]))
  U <- sapply(seq_len(p), function(j) margins_X[[j]](x[, j]))
  c_YX <- fit_copula(cbind(V, U), copula_method, weights, ...)
  c_Y <- fit_copula(V, copula_method, weights, ...)

  if (length(weights) == 0) {
    weights <- rep(1, n)
  }
  weights <- weights / mean(weights)
  w <- function(x) {
    u <- sapply(seq_len(p), function(j) margins_X[[j]](x[j]))
    u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
    Vu <- cbind(V, matrix(rep(u, each = n), n, p))
    c_YX(Vu) / c_Y(V) * weights
  }
  structure(
    list(
      copula_method = copula_method,
      margin_method = margin_method,
      dots = list(...),
      y = y,
      weights = weights,
      n = n,
      q = q,
      p = p,
      w = w
    ),
    class = "eecop"
  )
}

fit_margin <- function(x, method, weights) {
  switch(method,
    "normal" = fit_margin_normal(x, weights),
    "kde" = fit_margin_kde(x, weights)
  )
}

fit_copula <- function(u, method, weights, ...) {
  if (NCOL(u) == 1) {
    return(function(u) rep(1, NROW(u)))
  }
  switch(method,
    "vine" = fit_copula_vine(u, weights, ...),
    "normal" = fit_copula_normal(u, weights),
    "kde" = fit_copula_kde(u, weights),
    "bernstein" = fit_copula_bernstein(u, weights)
  )
}

get_psi <- function(type, y) {
  switch(type,
    "expectile" = get_psi_expectile(y),
    "quantile" = get_psi_quantile(y)
  )
}

#' Title
#'
#' @param object
#' @param x
#' @param type
#' @param t
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom assertthat is.scalar is.string
predict.eecop <- function(object, x, type = "expectile", t = 0.5, ...) {
  if ((NCOL(x) == 1) & (NROW(x) == object$p)) {
    x <- t(x)
  }
  x <- as.matrix(x)
  assert_that(
    ncol(x) == object$p,
    (is.string(type) & (type %in% c("expectile", "quantile"))) | is.function(type),
    is.numeric(t)
  )

  tol <- sd(object$y) / NROW(object$y)
  out <- apply(
    x, 1, predict_one_x,
    psi = get_psi(type, object$y),
    t = t,
    w = object$w,
    range = range(object$y),
    tol = tol
  )
  matrix(unlist(out), NROW(x), length(t), byrow = TRUE)
}

predict_one_x <- function(x, psi, t, w, range, tol) {
  w_x <- w(t(x))
  range <- range + c(-0.25, 0.25) * diff(range)
  lapply(t, predict_one_t, psi = psi, w_x = w_x, range = range, tol = tol)
}

predict_one_t <- function(t, psi, w_x, range, tol) {
  Eg <- function(theta) mean(psi(theta, t) * w_x)
  root_solve(Eg, range, tol)
}
