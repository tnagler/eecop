#' Bootstrapping eecop models
#'
#' Given a fitted model for the weight function, [eecop_boot()] generates a
#' number of bootstrap replicates from this model. A multiplier bootstrap
#' procedure is used. The result can be passed to [predict.eecop_boot()] to
#' generate bootstrapped predictions, or [conf_int()] to compute confidence
#' intervals directly.
#'
#' @param object a fitted [eecop] object.
#' @param n_boot number of bootstrap replicates.
#' @param rxi a function generating bootstrap multipliers; the function needs to
#'   take the number of samples as its first argument and return a numeric
#'   vector of this length. Default is [rexp()], which corresponds to the
#'   'Bayesian bootstrap'. Note that, for the theory to work,
#'   one should normalize so that `E(rxi) = var(rxi) = 1`, e.g.
#'   with `function(n) {xi <- rxi(n); return((xi - mean(xi)) / sd(xi) + 1)}`.
#' @param cores an integer for the number of cores to use; if `cores > 1`,
#' estimation will be parallelized
#' within over `n_boot` (using [parallel::parLapply()]).
#' @param x covariate values to predict on; must match the format used for
#'   fitting the `eecop()` model.
#' @param type either `"quantile"`, `"expectile"`, `"mean"`, or `"variance"`.
#' @param t a vector of quantile/expectile levels.
#' @param trafo a function with signature `function(y)` with `y` the
#'   response (vector or matrix). The function should return a
#'   vector of length `NROW(object$y)` or a matrix with `NROW(object$y)` rows.
#' @param conf confidence level.
#' @param ... unused.
#'
#' @return An objecvt of class `eccop_boot` containing the original [eecop]
#' object and bootstrap replicates.
#' @seealso [eecop()], [predict.eecop_boot()]
#' @export
#' @importFrom assertthat is.count
#' @importFrom parallel makeCluster stopCluster parLapply
#' @examples
#' # model with continuous variables
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#'
#' fit <- eecop(y, x)
#'
#' bs_fits <- eecop_boot(fit, n_boot = 2)
#' preds <- predict(bs_fits, x[1:3, ])
#' CI <- conf_int(bs_fits, x[1:3, ], type = "quantile", t = c(0.5, 0.9))
eecop_boot <- function(object, n_boot = 100, rxi = stats::rexp, cores = 1) {
  assert_that(inherits(object, "eecop"), is.count(n_boot))
  assert_that(is.function(rxi), length(rxi(5)) == 5, is.numeric(rxi(5)))
  assert_that(is.count(cores))

  if (cores > 1) {
    cl <- makeCluster(cores)
    on.exit(try(stopCluster(cl), silent = TRUE))
    lapply <- function(...) parLapply(cl, ...)
  }

  n <- object$n
  args <- list(
    y = object$y,
    x = object$x,
    copula_method = object$copula_method,
    margin_method = object$margin_method,
    weights = object$weights
  )
  args <- utils::modifyList(args, object$dots)
  do_one <- function(b) {
    xi <- rxi(n)
    xi <- xi / mean(xi)
    args$weights <- args$weights * xi
    do.call(eecop, args)
  }
  reps <- lapply(seq_len(n_boot), do_one)
  out <- list(orig = object, boot = reps)
  class(out) <- "eecop_boot"
  out
}


#' @rdname eecop_boot
#' @export
predict.eecop_boot <- function(object, x, type = "expectile", t = 0.5,
                               trafo = function(y) y, cores = 1, ...) {
  assert_that(inherits(object, "eecop_boot"))
  assert_that(is.count(cores))

  orig <- predict(object$orig, x = x, type = type, t = t, trafo = trafo)

  if (cores > 1) {
    cl <- makeCluster(cores)
    on.exit(try(stopCluster(cl), silent = TRUE))
    lapply <- function(...) parLapply(cl, ...)
  }

  boot <- lapply(
    object$boot,
    function(o) predict(o, x, type = type, t = t, trafo = trafo)
  )
  list(orig = orig, boot = boot)
}

#' @rdname eecop_boot
#' @export
conf_int <- function(object, x, type = "expectile", t = 0.5,
                     trafo = function(y) y,
                     conf = 0.9, cores = 1, ...) {
  assert_that(inherits(object, "eecop_boot"))
  preds <- predict(object,
    x = x, type = type, t = t,
    trafo = trafo, cores = cores
  )
  bdim <- dim(preds$boot[[1]])
  if (is.null(bdim)) bdim <- length(preds$boot[[1]])
  boot_arr <- array(unlist(preds$boot), dim = c(bdim, length(preds$boot)))

  alph <- (1 - conf) / 2
  fix <- seq_along(dim(boot_arr)[-1])
  low <- apply(boot_arr, fix, stats::quantile, probs = alph)
  up <- apply(boot_arr, fix, stats::quantile, probs = 1 - alph)
  mid <- apply(boot_arr, fix, mean)

  list(
    lower = 2 * preds$orig - up,
    estimate = preds$orig,
    upper =  2 * preds$orig - low
  )
}
