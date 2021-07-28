#' Prediction of quantiles or expectiles
#'
#' Predicts quantiles, expectiles, mean vectors, or covariate matrices from an
#' `eecop()` model by solving a weighted estimating equations as in Nagler and
#' Vatter (2020).
#'
#' @param object an `eecop` object.
#' @param x covariate values to predict on; must match the format used for
#'   fitting the `eecop()` model.
#' @param type either `"quantile"`, `"expectile"`, `"mean"`, or `"variance"`.
#' @param t a vector of quantile/expectile levels.
#' @param trafo a function with signature `function(y)` with `y` the
#'   response (vector or matrix). The function should return a
#'   vector of length `NROW(object$y)` or a matrix with `NROW(object$y)` rows.`
#' @param ... unused.
#'
#' @return For `"quantile"` and `"expectile"`: a matrix of predictions, each
#'   column corresponding to one `t` (in the order they were supplied to
#'   `predict()`. For `"mean"`, a matrix containing the predicted mean vectors
#'   in the rows. For  `"variance"` an 3-dimensional array containing the
#'   predicted variance/covariance matrices.
#' @export
#'
#' @references Nagler, T. and Vatter, T. (2020). Solving estimating equations
#'   with copulas. arXiv:1801.10576
#'
#' @examples
#' # univariate response
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#'
#' fit <- eecop(y, x)
#' predict(fit, x[1:3, ], t = c(0.5, 0.9), type = "quantile")
#' predict(fit, x[1:3, ], t = c(0.5, 0.9), type = "expectile")
#'
#' # multivariate response
#' y <- cbind(y1 = y, y2 = y + rnorm(100))
#' fit <- eecop(y, x)
#'
#' predict(fit, x[1:3, ], type = "mean")
#' predict(fit, x[1:3, ], type = "variance")
#' predict(fit, x[1:3, ], type = "quantile", trafo = function(y) y[, 1] + y[, 2])
#' @importFrom assertthat is.scalar is.string
#' @importFrom stats predict
predict.eecop <- function(object,
                          x,
                          type = "expectile",
                          t = 0.5,
                          trafo = function(y) y,
                          ...) {
  assert_that(is.string(type))
  y <- as.data.frame(trafo(object$y))
  x <- process_x_new(object, x)
  switch(type,
    "mean" = predict_mean(object, y, x),
    "variance" = predict_variance(object, y, x),
    "expectile" = predict_expectile(object, y, x, t),
    "quantile" = predict_quantile(object, y, x, t),
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

predict_mean <- function(object, y, x) {
  pred <- vapply(
    seq_len(nrow(x)),
    function(i) {
      w <- object$w(x[i, , drop = FALSE])
      colSums(y * w / sum(w))
    },
    numeric(ncol(y))
  )
  if (ncol(y) > 1) {
    pred <- t(pred)
  }
  pred
}

predict_variance <- function(object, y, x) {
  vapply(
    seq_len(nrow(x)),
    function(i) {
      w <- object$w(x[i, , drop = FALSE])
      stats::cov.wt(y, wt = w / sum(w))$cov
    },
    matrix(1, ncol(y), ncol(y))
  )
}

predict_uniroot <- function(object, y, x, t, idfun) {
  out <- sapply(
    seq_len(nrow(x)),
    function(i, ...) predict_uni_x(x[i, , drop = FALSE], ...),
    psi = idfun,
    t = t,
    w = object$w,
    range = range(y),
    tol = sd(y) / object$n
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

predict_quantile <- function(object, y, x, t = 0.5) {
  if (ncol(y) > 1) {
    stop("can't predict quantiles for multivariate response.")
  }
  y <- y[[1]]
  assert_that(is.numeric(t), all((0 < t) & (t < 1)))

  idfun <- function(theta, t) {
    (y <= theta) - t
  }
  predict_uniroot(object, y, x, t, idfun)
}

predict_expectile <- function(object, y, x, t = 0.5) {
  if (ncol(y) > 1) {
    stop("can't predict expectiles for multivariate response.")
  }
  y <- y[[1]]
  assert_that(is.numeric(t), all((0 < t) & (t < 1)))

  idfun <- function(theta, t) {
    t * (y - theta) * (y >= theta) - (1 - t) * (theta - y) * (y < theta)
  }
  predict_uniroot(object, y, x, t, idfun)
}



#' Generic solver for copula-based estimating equations
#'
#' Solves an estimating equation based on a fitted [eecop] model and
#' user-supplied identifying function.
#'
#' @param object a fitted [eecop] object.
#' @param x covariate values to predict on; must match the format used for
#'   fitting the `eecop()` model.
#' @param idfun a function with signature `function(y, theta)` with `y` the
#'   response (vector or matrix) and `theta` the parameter of interest. The
#'   function should return a vector of length `NROW(object$y)` or a matrix with
#'   `NROW(object$y)` rows.
#' @param theta_start starting values for optimizing the parameter of interest.
#' @param ... further arguments passed to [optim()].
#'
#' @return The optimal parameter `theta`.
#' @export
#'
#' @examples
#' ## fit dummy model
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#' fit <- eecop(y, x)
#'
#' ## identifying function for 0.5 and 0.9 quantiles
#' idfun <- function(y, theta) {
#'   t <- c(0.5, 0.9)
#'   gmat <- matrix(NA, NROW(y), 2)
#'   for (j in 1:2) gmat[, j] <- (y <= theta[j]) - t[j]
#'   gmat
#' }
#'
#' ## solve estimating equation
#' solve_eecop(fit, x[1:3, ], idfun = idfun, theta_start = rep(0, 2))
solve_eecop <- function(object, x, idfun, theta_start, ...) {
  lapply(
    seq_len(nrow(x)),
    function(i, ...) {
      w <- object$w(x[i, , drop = FALSE]) * object$weights
      w <- w / sum(w)
      Eg2 <- function(theta) sqrt(sum(colSums(idfun(object$y, theta) * w)^2))
      stats::optim(theta_start, fn = Eg2, ...)$par
    },
    ...
  )
}
