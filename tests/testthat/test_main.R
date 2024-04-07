
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
    wgcm1 <- wgcm(Y, X, Z)
  })
})

test_that("comet", {
  set.seed(42)
  data("mtcars")
  expect_no_error({
    comet(mpg ~ cyl | disp, data = mtcars)
    comet(factor(vs) ~ cyl | disp, data = mtcars)
    comet(mpg ~ cyl | disp, data = mtcars, test = "pcm")
    comet(factor(vs) ~ cyl | disp, data = mtcars, test = "pcm")
    comet(mpg ~ cyl | disp, data = mtcars, test = "wgcm")
    comet(factor(vs) ~ cyl | disp, data = mtcars, test = "wgcm")
  })
  expect_error({
    comet(factor(mpg) ~ cyl | disp, data = mtcars)
    comet(factor(mpg) ~ cyl | disp, data = mtcars, test = "pcm")
    comet(factor(mpg) ~ cyl | disp, data = mtcars, test = "wgcm")
  })
})
