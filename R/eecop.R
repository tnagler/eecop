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
    copula_method %in% c("vine", "normal", "kde", "trafo_kde"),
    margin_method %in% c("kde", "normal"),
    is.numeric(weights)
  )

  if (any(sapply(y, function(yy) is.factor(y) & !is.ordered(y)))) {
    stop("factor-valued response not allowed.")
  }
  if (copula_method == "kde") {
    margin_method <- "void"
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
                     method = copula_method,
                     weights = weights,
                     var_types = c(var_types_Y, var_types_X),
                     ...
  )
  c_Y <- fit_copula(V, copula_method, weights, var_types = var_types_Y, ...)

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
         "kde" = fit_margin_kde(x, weights),
         "void" = fit_margin_void(x, weights)
  )
}

fit_copula <- function(u, method, weights, var_types, ...) {
  if (length(var_types) == 1) {
    return(function(u) rep(1, NROW(u)))
  }
  switch(method,
         "vine" = fit_copula_vine(u, weights, ...),
         "normal" = fit_copula_normal(u, weights),
         "kde" = fit_copula_kde(u, weights, ...),
         "trafo_kde" = fit_copula_trafo_kde(u, weights, ...)
  )
}

get_psi <- function(type, y) {
  switch(type,
         "expectile" = get_psi_expectile(y),
         "quantile" = get_psi_quantile(y)
  )
}

#' Prediction of quantiles or expectiles
#'
#' Predicts quantiles or expectiles from an `eecop()` model by solving
#' a weighted estimating equations as in Nagler and Vatter (2020).
#'
#' @param object an `eecop` object.
#' @param x covariate values to predict on; must match the format used for
#' fitting the `eecop()` model.
#' @param type either `"quantile"` or `"expectile"`.
#' @param t a vector of quantile/expectile levels.
#' @param ... unused.
#'
#' @return
#' A matrix of predictions, each column corresponding to one `t` (in the order
#' they were supplied to `predict()`.
#' @export
#'
#' @references
#' Nagler, T. and Vatter, T. (2020). Solving estimating equations with copulas.
#' arXiv:1801.10576
#'
#' @examples
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#'
#' fit <- eecop(y, x)
#' predict(fit, x, t = c(0.5, 0.9), type = "quantile")
#' predict(fit, x, t = c(0.5, 0.9), type = "expectile")
#' @importFrom assertthat is.scalar is.string
predict.eecop <- function(object, x, type = "expectile", t = 0.5) {
  assert_that(is.string(type))
  switch(
    type,
    "mean" = predict_mean(object, x),
    "variance" = predict_variance(object, x),
    "expectile" = predict_expectile(object, x, t = t),
    "quantile" = predict_quantile(object, x, t = t),
    stop(paste0("type '", type, "' not implemented"))
  )
}

process_x_new <- function(object, x) {
  if ((NCOL(x) == 1) & (NROW(x) == object$p)) {
    x <- t(x)
  }
  x <- as.data.frame(x)
  assert_that(ncol(x) == object$p)
  x
}

predict_mean <- function(object, x) {
  x <- process_x_new(object, x)
  pred <- vapply(
    seq_len(nrow(x)),
    function(i) {
      w <- object$w(x[i, , drop = FALSE]) * object$weights
      colMeans(object$y * w / sum(w))
    },
    numeric(object$q)
  )
  if (object$q > 1) {
    pred <- t(pred)
  }
  pred
}


predict_variance <- function(object, x) {
  x <- process_x_new(object, x)
  vapply(
    seq_len(nrow(x)),
    function(i) {
      w <- object$w(x[i, , drop = FALSE]) * object$weights
      cov.wt(object$y, wt = w / sum(w))$cov
    },
    matrix(1, object$q, object$q)
  )
}

predict_uniroot <- function(object, x, t, idfun) {
  x <- process_x_new(object, x)
  tol <- sd(object$y[[1]]) / object$n
  out <- sapply(
    seq_len(nrow(x)),
    function(i, ...) predict_uni_x(x[i, , drop = FALSE], ...),
    psi = idfun,
    t = t,
    w = object$w,
    range = range(object$y),
    tol = sd(object$y[[1]]) / object$n
  )
  matrix(unlist(out), NROW(x), length(t), byrow = TRUE)
}


predict_uni_x <- function(x, psi, t, w, range, tol) {
  w_x <- w(x)
  w_sel <- which(!is.nan(w_x))
  if (!length(w_sel)) {
    return(lapply(t, function(tt) NA))
  }
  range <- range + c(-0.25, 0.25) * diff(range)
  lapply(t, predict_uni_t,
         psi = psi, w_x = w_x[w_sel],
         w_sel = w_sel, range = range, tol = tol
  )
}

predict_uni_t <- function(t, psi, w_x, w_sel, range, tol) {
  Eg <- function(theta) mean(psi(theta, t)[w_sel] * w_x)
  root_solve(Eg, range, tol)
}

predict_quantile <- function(object, x, t = 0.5) {
  if (object$q > 1)
    stop("can't predict quantiles for multivariate response.")
  x <- process_x_new(object, x)
  assert_that(is.numeric(t), all((0 < t) & (t < 1)))

  idfun <- function(theta, t) {
    (object$y[[1]] <= theta) - t
  }
  predict_uniroot(object, x, t, idfun)
}

predict_expectile <- function(object, x, t = 0.5) {
  if (object$q > 1)
    stop("can't predict expectiles for multivariate response.")
  x <- process_x_new(object, x)
  assert_that(is.numeric(t), all((0 < t) & (t < 1)))

  y <- object$y[[1]]
  idfun <-  function(theta, t) {
    t * (y - theta) * (y >= theta) - (1 - t) * (theta - y) * (y < theta)
  }
  predict_uniroot(object, x, t, idfun)
}

predict_gen <- function(object, x, idfun, theta_start, ...) {
  x <- process_x_new(object, x)
  lapply(
    seq_len(nrow(x)),
    function(i, ...) {
      browser()
      w <- object$w(x[i, , drop = FALSE]) * object$weights
      w <- w / sum(w)
      Eg2 <- function(theta) sqrt(sum(colSums(idfun(object$y, theta) * w)^2))
      optim(theta_start, fn = Eg2, ...)$par
    },
    ...
  )
}


