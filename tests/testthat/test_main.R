
test_that("gcm and pcm work", {
  expect_no_error({
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    gcm1 <- gcm(Y, X, Z)
    pcm1 <- pcm(Y, X, Z)
  })
})

test_that("comet", {
  set.seed(42)
  data("mtcars")
  expect_no_error(comet(mpg ~ cyl | disp, data = mtcars))
  expect_no_error(comet(factor(vs) ~ cyl | disp, data = mtcars))
  expect_error(comet(factor(mpg) ~ cyl | disp, data = mtcars))
})
