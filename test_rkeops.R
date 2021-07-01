library(rkeops)
library(ks)

# to run computation on CPU (default mode)
use_cpu()
# to run computations on GPU (to be used only if relevant)
use_gpu()

# Pre-compile the KeOps reductions
precompile_reductions <- function(d_max = 10) {
  for (d in 1:d_max) {

    # KeOps formula
    formula <- "Sum_Reduction(Exp(-SqDist(x, y)), 1)"
    # KeOps inputs
    args <- c(
      sprintf("x = Vi(%d)", d),
      sprintf("y = Vj(%d)", d)
    )

    # KeOps kernel compilation
    reduction <- rkeops::keops_kernel(formula, args)
  }
  invisible(d_max)
}

# KDE implementation with KeOps
kde_keops <- function(x) {

  # Dimensions
  d <- ncol(x)
  n <- nrow(x)

  # Bandwidth
  bw <- n^(-2 / (4 + d)) * cov(x)
  bw_ir <- solve(chol(bw)) / sqrt(2)

  # Rotate the input
  z <- x %*% bw_ir

  # KeOps formula
  formula <- "Sum_Reduction(Exp(-SqDist(x, y)), 1)"
  # KeOps inputs
  args <- c(
    sprintf("x = Vi(%d)", d),
    sprintf("y = Vj(%d)", d)
  )

  # KeOps kernel compilation
  reduction <- rkeops::keops_kernel(formula, args)

  # Compute the reduction
  res <- as.numeric(reduction(list(z, z)))

  # Scale and return result
  return(res / (n * sqrt(2^d * pi^d * det(bw))))
}

# KDE implementation with KS
kde_ks <- function(x) {

  # Dimensions
  d <- ncol(x)
  n <- nrow(x)

  # Bandwidth
  bw <- n^(-2 / (4 + d)) * cov(x)
  bws <- replicate(n, bw, simplify = FALSE)
  bws <- do.call(rbind, bws)

  # Compute the kde and return
  weights <- rep(1 / n, n)
  if (d == 1) {
    res <- ks::dnorm.mixt(
      x = x, mus = x,
      sigmas = rep(sqrt(bws), n),
      props = weights
    )
  } else {
    res <- ks::dmvnorm.mixt(
      x = x, mus = x,
      Sigmas = bws,
      props = weights
    )
  }
  return(res)
}

# This takes around d_max * 15 seconds
precompile_reductions()

# Some 1-d data
n <- 1000
d <- 1
x <- matrix(rnorm(n * d), nrow = n)

# Compare the compute times
system.time(res1 <- kde_keops(x))
system.time(res2 <- kde_ks(x))
res3 <- density(x) # Sanity check

# Everything loooks ok
xx <- sort(x, index = TRUE)
par(mfrow = c(1, 1))
plot(xx$x, dnorm(xx$x), type = "l", xlab = "x", ylab = "density")
lines(xx$x, res1[xx$ix], col = "blue")
lines(xx$x, res2[xx$ix], col = "red")
lines(res3$x, res3$y, col = "green")

# Some 3-d data
d <- 3
x <- matrix(rnorm(n * d), nrow = n)

# Compare the compute times
system.time(res1 <- kde_keops(x))
system.time(res2 <- kde_ks(x))

# Sanity check
par(mfrow = c(1, 3), pty = "s")
truth <- apply(dnorm(x), 1, prod)
for (j in 1:3) {
  xx <- sort(x[, j], index = TRUE)
  plot(xx$x, truth[xx$ix],
    xlab = "x", ylab = "density",
    main = sprintf("X%d", j)
  )
  points(xx$x, res1[xx$ix], col = "blue")
  points(xx$x, res2[xx$ix], col = "red")
}


# Just for fun
tt <- matrix(NA, 10, 2)
for (d in 1:10) {
  x <- matrix(rnorm(n * d), nrow = n)

  # Compare the compute times
  tt[d, 1] <- system.time(res1 <- kde_keops(x))[3]
  tt[d, 2] <- system.time(res2 <- kde_ks(x))[3]
}
matplot(1:10, tt,
  type = "l",
  xlab = "d", ylab = "time", main = "n = 1000"
)

d <- 10
nn <- 10^(1:5)
t <- vector(length = 5)
for (j in 1:5) {
  n <- nn[j]
  x <- matrix(rnorm(n * d), nrow = n)
  t[j] <- system.time(res <- kde_keops(x))[3]
}

n <- 10^4
d <- 10
x <- matrix(rnorm(n * d), nrow = n)

use_cpu(ncore = 1)
system.time(res <- kde_keops(x))
system.time(res <- kde_keops(x))

use_cpu(ncore = 0)
system.time(res <- kde_keops(x))
system.time(res <- kde_keops(x))

clean_rkeops()
use_gpu()
compile4gpu()
system.time(res <- kde_keops(x))
system.time(res <- kde_keops(x))