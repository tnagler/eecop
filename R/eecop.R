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
  y <- as.data.frame(y)
  x <- as.data.frame(x)
  assert_that(
    nrow(y) == nrow(x),
    is.string(copula_method),
    is.string(margin_method),
    copula_method %in% c("vine", "normal", "kde", "bernstein"),
    margin_method %in% c("kde", "normal"),
    is.numeric(weights)
  )

  if (any(sapply(y, function(v) is.factor(v) & !is.ordered(v))))
    stop("factor-valued response not allowed.")
  x <- rvinecopulib:::expand_factors(x)

  q <- ncol(y)
  p <- ncol(x)
  n <- nrow(x)

  var_types_Y <- ifelse(sapply(y, is.ordered), "d", "c")
  var_types_X <- ifelse(sapply(x, is.ordered), "d", "c")

  if ("d" %in% c(var_types_X, var_types_Y)) {
    if (margin_method == "normal")
      stop("normal margins can't be used with discrete data.")
    if (copula_method != "vine")
      stop('only copula_method = "vine" allowed with discrete data.')
  }

  margins_Y <- lapply(
    seq_len(q),
    function(j) fit_margin(y[, j], margin_method, weights)
  )
  margins_X <- lapply(
    seq_len(p),
    function(j) fit_margin(x[, j], margin_method, weights)
  )

  V <- compute_pseudo_obs(y, margins_Y)
  U <- compute_pseudo_obs(x, margins_X)
  c_YX <- fit_copula(combine_margins(V, U, q, p),
                     copula_method, weights,
                     var_types = c(var_types_Y, var_types_X),
                     ...)
  c_Y <- fit_copula(V, copula_method, weights,
                    var_types = var_types_Y,
                    ...)

  if (length(weights) == 0) {
    weights <- rep(1, n)
  }
  weights <- weights / mean(weights)
  w <- function(x) {
    u <- compute_pseudo_obs(x, margins_X)
    u <- matrix(rep(u, each = n), n, ncol(u))
    Vu <- combine_margins(V, u, q, p)
    c_YX(Vu) / c_Y(V) * weights
  }
  structure(
    list(
      copula_method = copula_method,
      margin_method = margin_method,
      dots = list(...),
      var_types_Y = var_types_Y,
      var_types_X = var_types_X,
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
         "kde" = fit_copula_kde(u, weights, ...),
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
  x <- as.data.frame(x)
  assert_that(
    ncol(x) == object$p,
    (is.string(type) & (type %in% c("expectile", "quantile"))) | is.function(type),
    is.numeric(t)
  )

  tol <- max(apply(object$y, 2, sd)) / NROW(object$y)
  out <- sapply(
    seq_len(nrow(x)),
    function(i, ...) predict_one_x(x[i, , drop = FALSE], ...),
    psi = get_psi(type, object$y),
    t = t,
    w = object$w,
    range = range(object$y),
    tol = tol
  )
  matrix(unlist(out), NROW(x), length(t), byrow = TRUE)
}

predict_one_x <- function(x, psi, t, w, range, tol) {
  w_x <- w(x)
  w_sel <- which(!is.nan(w_x))
  range <- range + c(-0.25, 0.25) * diff(range)
  lapply(t, predict_one_t,
         psi = psi, w_x = w_x[w_sel],
         w_sel = w_sel, range = range, tol = tol
  )
}

predict_one_t <- function(t, psi, w_x, w_sel, range, tol) {
  Eg <- function(theta) mean(psi(theta, t)[w_sel] * w_x)
  root_solve(Eg, range, tol)
}
