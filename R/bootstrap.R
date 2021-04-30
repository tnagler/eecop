#' Bootstrapping eecop models
#'
#' Given a fitted model for the weight function, [bootstrap.eecop] generates a
#' number of bootstrap replicates from this model. A multiplier bootstrap
#' procedure is used.
#'
#' @param object a fitted [eecop] object.
#' @param n_boot number of bootstrap replicates.
#' @param rxi a function generating bootstrap multipliers; the function needs to
#'   take the number of samples as its first argument and return a numeric
#'   vector of this length. Default is [rexp()], which corresponds to the
#'   'Bayesian bootstrap'.
#'
#' @return A list of [eecop] objects.
#' @seealso [eecop()], [predict.bootstrap_list()]
#' @export
#' @importFrom assertthat is.count
#' @examples
#' # model with continuous variables
#' x <- matrix(rnorm(200), 100, 2)
#' y <- rowSums(x) + rnorm(100)
#'
#' fit <- eecop(y, x)
#'
#' bs_fits <- bootstrap(fit, n_boot = 2)
#' preds <- predict(bs_fits)
bootstrap.eecop <- function(object, n_boot = 100, rxi = rexp) {
  assert_that(inherits(object, "eecop"), is.count(n_boot))
  assert_that(is.function(rxi), length(rxi(5)) == 5, is.numeric(rxi(5)))

  n <- object$n
  args <- list(
    y = object$y,
    x = object$x,
    copula_method = object$copula_method,
    margin_method = object$margin_method,
    weights = object$weights
  )
  args <- modifyList(args, object$dots)
  out <- lapply(seq_len(n_boot), function(b) {
    xi <- rxi(n)
    xi <- (xi - mean(xi)) / sd(xi) + 1
    args$weights <- args$weights * xi
    do.call(eecop, args)
  })
  class(out) <- "eecop_list"
  out
}
