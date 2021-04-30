#' Copula regression models based on estimating equations
#'
#' Implements the copula regression estimators of Nagler and Vatter (2020). A
#' model for marginal distributions and copula between response and covariates
#' is estimated. Predictions for quantiles or expectiles can then be derived
#' from solving a weighted estimating equations.
#'
#' @param y vector with numeric response values.
#' @param x vector, matrix, or data.frame with covariate values (rows are
#'   observations).
#' @param copula_method method for estimating the copula(s); one of `vine",
#'   "normal", "kde"` for vine copula, Gaussian copula, and transformation
#'   kernel density method, respectively.
#' @param margin_method method for estimating marginal distributions; one of
#'   `"kde", "normal"` for kernel density or Gaussian margins, respectively.
#' @param weights optional; a vector of weights for each observation.
#' @param ... further arguments passed to `rvinecopulib::vinecop()`.
#'
#' @details Both `y` and `x` may contain discrete variables, which must be
#' passed as `ordered()` or `factor()` variables.
#'
#' @return An object of class `eecop`. Use `predict.eecop()` to predict
#'   quantiles or expectiles. For other estimating equations, the weights
#'   \eqn{w(x) = c_{YX}(Y, x)/ c_Y(Y)} can be computed from `object$w(x)`.
#'
#' @seealso `predict.eecop()`
#' @export
#'
#' @importFrom assertthat assert_that is.string
#'
#' @references Nagler, T. and Vatter, T. (2020). Solving estimating equations
#' with copulas. arXiv:1801.10576
#'
#' @examples
#' # model with continuous variables
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#'
#' fit <- eecop(y, x)
#'
#' predict(fit, x, t = c(0.5, 0.9), type = "quantile")
#' predict(fit, x, t = c(0.5, 0.9), type = "expectile")
#'
#' # model with discrete covariates
#' x <- as.data.frame(matrix(rbinom(200, 5, 0.3), 100, 2))
#' y <- rowSums(x) + rnorm(100)
#' for (k in 1:2) {
#'   x[, k] <- ordered(x[, k], levels = 0:5)
#' }
#'
#' fit <- eecop(y, x)
#'
#' predict(fit, x, t = c(0.5, 0.9), type = "quantile")
#' predict(fit, x, t = c(0.5, 0.9), type = "expectile")
eecop <- function(y, x, copula_method = "vine", margin_method = "kde",
                  weights = numeric(), ...) {
  y <- as.data.frame(y)
  x <- as.data.frame(x)
  assert_that(
    nrow(y) == nrow(x),
    is.string(copula_method),
    is.string(margin_method),
    copula_method %in% c("vine", "normal", "kde"),
    margin_method %in% c("kde", "normal"),
    is.numeric(weights)
  )

  if (any(sapply(y, function(yy) is.factor(y) & !is.ordered(y)))) {
    stop("factor-valued response not allowed.")
  }


  x <- rvinecopulib:::expand_factors(x)

  q <- ncol(y)
  p <- ncol(x)
  n <- nrow(x)

  var_types_Y <- ifelse(sapply(y, is.ordered), "d", "c")
  var_types_X <- ifelse(sapply(x, is.ordered), "d", "c")

  if ("d" %in% c(var_types_Y, var_types_X)) {
    if (margin_method == "normal") {
      stop("normal margins can't be used with discrete data.")
    }
    if (copula_method != "vine") {
      stop('only copula_method = "vine" allowed with discrete data.')
    }
  }

  if (copula_method == "kde") {
    margin_method <- "void"
  }
  margins_Y <- lapply(
    seq_len(q),
    function(j) fit_margin(y[, j], margin_method, weights)
  )
  margins_X <- lapply(
    seq_len(p),
    function(j) fit_margin(x[, j], margin_method, weights)
  )
  if (margin_method == "void") {
    margin_method <- "kde"
  }

  V <- compute_pseudo_obs(y, margins_Y)
  U <- compute_pseudo_obs(x, margins_X)

  w_model <- fit_w(V, U,
    method = copula_method,
    weights = weights,
    var_types = c(var_types_Y, var_types_X),
    ...
  )

  if (length(weights) == 0) {
    weights <- rep(1, n)
  }
  weights <- weights / mean(weights)
  w <- function(x) {
    u <- compute_pseudo_obs(x, margins_X)
    u <- matrix(rep(u, each = n), n, ncol(u))
    w_model(u) * weights
  }
  structure(
    list(
      copula_method = copula_method,
      margin_method = margin_method,
      dots = list(...),
      var_types_Y = var_types_Y,
      var_types_X = var_types_X,
      y = y,
      x = x,
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
    "kde" = fit_margin_kde(x, weights),
    "void" = fit_margin_void(x, weights)
  )
}

fit_w <- function(v, u, method, weights, var_types, ...) {
  if (length(var_types) == 1) {
    return(function(u) rep(1, NROW(u)))
  }
  switch(method,
    "vine" = fit_w_vine(v, u, weights, ...),
    "normal" = fit_w_normal(v, u, weights),
    "kde" = fit_w_kde(v, u, weights, ...)
  )
}
