#' @noRd
#' @importFrom stats uniroot
root_solve <- function(f, range, tol = .Machine$double.eps^0.25) {
  root <- tryCatch(uniroot(f, range, tol = tol)$root, error = function(e) e)

  if (any(class(root) == "error")) {
    x0 <- mean(range) - (mean(range) - min(range)) / 20
    x1 <- mean(range) + (max(range) - mean(range)) / 20
    root <- root_solve_secant(f, x0, x1, tol)
  }
  return(root)
}

root_solve_secant <- function(fun, x0, x1,
                              tol = .Machine$double.eps^0.25,
                              maxiter = 1000) {
  i <- 1
  for (i in 1:maxiter) {
    x2 <- x1 - fun(x1) * (x1 - x0) / (fun(x1) - fun(x0))
    if (abs(fun(x2)) < tol) {
      return(x2)
    }
    x0 <- x1
    x1 <- x2
  }
  warning("max number of iterations reached")
  return(x2)
}
