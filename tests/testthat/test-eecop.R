set.seed(5)

n <- 50
d <- 3

x <- matrix(rnorm(n * d), n, d)
y <- rowSums(x) + rnorm(n)

xx <- lapply(seq_len(ncol(x)), function(j) as.ordered(rbinom(n, 4, 0.2)))
xx <- as.data.frame(xx)

w <- rexp(n)

test_that("continuous fitting works", {
  expect_silent(eecop(y, x, margin_method = "normal"))
  expect_silent(eecop(y, x, copula_method = "kde"))
  expect_silent(eecop(y, x, copula_method = "normal"))
  expect_silent(eecop(y, x, mult = 2))
  expect_silent(eecop(y, x, weights = rexp(length(y))))
  expect_error(eecop(y, x, weights = rexp(2)))
  expect_error(eecop(y, x, copula_method = "bernstein"))
  expect_error(eecop(c(y, y), x, copula_method = "bernstein"))
})

test_that("discrete fitting works", {
  expect_silent(eecop(y, xx))
  expect_silent(eecop(y, xx, weights = w))
  expect_error(eecop(y, xx, margin_method = "kde", copula_method = "kde"))
  y <- as.ordered(y)
  expect_silent(eecop(y, xx))
})


test_that("prediction works", {
  for (copula_method in c("vine", "kde", "normal")) {
    fit <- eecop(y, x, copula_method = copula_method)
    expect_equal(dim(predict(fit, x, t = 1:3 / 4)), c(nrow(x), 3))
    expect_equal(dim(predict(fit, x, "quantile", t = 1:3 / 4)), c(nrow(x), 3))
  }

  fit <- eecop(y, xx)
  expect_equal(dim(predict(fit, xx, t = 1:3 / 4)), c(nrow(x), 3))
  expect_equal(dim(predict(fit, xx, "quantile", t = 1:3 / 4)), c(nrow(x), 3))

  y <- cbind(y, y + rnorm(length(y)))
  fit <- eecop(y, x)
  expect_equal(dim(predict(fit, x, "mean")), c(nrow(x), 2))
  expect_equal(dim(predict(fit, x, "variance", t = 1:3 / 4)), c(2, 2, nrow(x)))
})


test_that("bootstrap works", {
  fit <- eecop(y, x)
  fit_bs <- bootstrap(fit, n_boot = 2)
  expect_length(fit_bs, 2)
  expect_true(all(sapply(fit_bs, class) == "eecop"))

  pred_exp <- predict(fit_bs, x, t = 1:3 / 4)
  expect_equal(dim(pred_exp[[1]]), c(nrow(x), 3))
  expect_length(pred_exp, 2)
})
