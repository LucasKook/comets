
test_that(".ranger works", {
  set.seed(12)
  tn <- 3e2
  dat <- data.frame(bin = factor(sample(0:1, tn, TRUE)),
                    ord = ordered(sample(7:10, tn, TRUE)),
                    mcc = factor(sample(11:15, tn, TRUE)),
                    num = rnorm(tn),
                    x = rnorm(tn))
  lapply(colnames(dat)[-ncol(dat)], \(resp) {
    fm <- reformulate("x", resp)
    rf <- .ranger(fm, data = dat)
    rr <- residuals.ranger(rf)
    expect_lt(abs(mean(rr)), 0.05)
  })
})

test_that("gcm pcm data types", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    X <- rnorm(tn)
    Z <- rnorm(tn)
    Y <- rnorm(tn)
    gcm1 <- gcm(Y, X, Z)
    pcm1 <- pcm(Y, X, Z)
    gcm2 <- gcm(Y, cbind(X, rnorm(tn)), cbind(Z, rnorm(tn)))
    gcm2 <- gcm(Y, as.data.frame(X), as.data.frame(Z))
  })
})

test_that("check data works", {
  tmp <- rnorm(10)
  expect_equal(NCOL(cY <- .check_data(tmp, "Y")), 1)
  expect_null(colnames(cY))
  expect_equal(ncol(cX <- .check_data(tmp, "X")), 1)
  expect_equal(colnames(cX), "X1")
  expect_equal(NCOL(.check_data(cbind(tmp, tmp), "Z")), 2)
})

test_that("wGCM with different regressions", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    wgcm1 <- wgcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "lasso")
    wgcm2 <- wgcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "rf")
    wgcm3 <- wgcm(Y, X, Z, reg_XonZ = "ridge", reg_YonZ = "ridge")
    wgcm4 <- wgcm(Y, X, Z, reg_XonZ = "qrf", reg_YonZ = "qrf")
    tmp <- plot(wgcm4, plot = FALSE)
  })
})

test_that("GCM with different regressions", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    gcm1 <- gcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "lasso")
    gcm2 <- gcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "rf",
                args_YonZ = list(mtry = 2))
    gcm3 <- gcm(Y, X, Z, reg_XonZ = "ridge", reg_YonZ = "ridge")
    gcm4 <- gcm(Y, X, Z, reg_XonZ = "qrf", reg_YonZ = "qrf")
    gcm5 <- gcm(Y, X, Z, reg_XonZ = "postlasso", reg_YonZ = "postlasso")
    gcm6 <- gcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm")
    tmp <- plot(gcm6, plot = FALSE)
  })
})

test_that("PCM with different regressions", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    pcm1 <- pcm(Y, X, Z, reg_YonXZ = "lasso", reg_YonZ = "lasso")
    pcm2 <- pcm(Y, X, Z, reg_YonXZ = "rf", reg_YonZ = "lasso")
    pcm3 <- pcm(Y, X, Z, reg_YonXZ = "ridge", reg_YonZ = "ridge")
    pcm4 <- pcm(Y, X, Z, reg_YonXZ = "postlasso", reg_YonZ = "postlasso")
    pcm5 <- pcm(Y, X, Z, reg_YonXZ = "lrm", reg_YonZ = "lrm")
    tmp <- plot(pcm5, plot = FALSE)
  })
})

test_that("Multi-dimensional GCM works", {
  set.seed(12)
  tn <- 3e2
  set.seed(12)
  X <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(X) <- c("X1", "X2")
  Z <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- cbind(rnorm(tn), rnorm(tn), rnorm(tn))
  expect_no_error(gcm1 <- gcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "lasso"))
  tmp <- plot(gcm1, plot = FALSE)
  ### With multi-level factor as Y
  expect_no_error(comet(Species ~ Sepal.Length | Sepal.Width, data = iris))
  expect_error(comet(Species ~ Sepal.Length | Sepal.Width, data = iris, test = "pcm"))
  expect_error(comet(Species ~ Sepal.Length | Sepal.Width, data = iris, test = "wgcm"))
})

test_that("TRAM GCM works with coxph and survforest", {
  library("survival")
  data("GBSG2", package = "TH.data")
  y <- Surv(GBSG2$time, GBSG2$cens)
  x <- model.matrix(~ 0 + horTh, data = GBSG2)[, 2]
  m <- cox(y, x)
  z <- model.matrix(~ 0 + age, data = GBSG2)
  expect_no_error(tgcm <- gcm(y, x, z, reg_YonZ = "cox"))
  expect_no_error(tgcm <- gcm(y, x, z, reg_YonZ = "survforest"))
  expect_no_error(comet(Surv(time, cens) ~ horTh | age, data = GBSG2,
                        reg_YonZ = "cox"))
})

test_that("coin for tests", {
  expect_no_error({
    library("coin")
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- cbind(rowSums(X) + rnorm(tn), rowSums(X) + rnorm(tn))
    gcm(Y, X, Z, type = "max", coin = TRUE, cointrol = list(
      distribution = approximate(99)))
    gcm(Y, X, Z, type = "quadratic", coin = TRUE, cointrol = list(
      distribution = "asymptotic"))
    gcm(Y[, 1], X[, 1], Z, type = "scalar", coin = TRUE, cointrol = list(
      distribution = "asymptotic"))
  })
})

test_that("equivalence test on different scales", {
  n <- 150
  X <- rnorm(n)
  Z <- matrix(rnorm(2 * n), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- X^2 + Z[, 2] + rnorm(n)
  expect_no_error(plm_equiv_test(Y, X, Z, from = -1, to = 1))
  expect_no_error(plm_equiv_test(Y, X, Z, from = -1, to = 1, scale = "cov"))
  expect_no_error(plm_equiv_test(Y, X, Z, from = -1, to = 1, scale = "cor"))
  expect_error(plm_equiv_test(Y, cbind(X, X), Z, from = -1, to = 1, scale = "cor"))
})

test_that("PCM works w/ and w/o replacement and fixed indices", {
  n <- 150
  X <- rnorm(n)
  Z <- matrix(rnorm(2 * n), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- X^2 + Z[, 2] + rnorm(n)
  expect_no_error(t1 <- pcm(Y, X, Z, indices = 1:75))
  expect_equal(t1$check.data$id, 76:n)
  expect_error(pcm(Y, X, Z, indices = 1:75, rep = 2))
})

test_that("fitted models can be returned", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- cbind(rowSums(X) + rnorm(tn), rowSums(X) + rnorm(tn))
    gcm <- gcm(Y, X, Z, return_fitted_models = TRUE)
    expect_length(tmp <- gcm$models[["reg_YonZ"]], NCOL(Y))
    expect_s3_class(tmp[[1]], c("rf", "ranger"))
    expect_length(gcm$models[["reg_XonZ"]], NCOL(X))
    pcm <- pcm(Y[, 1], X, Z, return_fitted_models = TRUE)
    expect_s3_class(pcm$models$reg_YonXZ, c("rf", "ranger"))
    wgcm <- wgcm(Y[, 1], X, Z, return_fitted_models = TRUE)
    expect_s3_class(wgcm$models[[1]], c("rf", "ranger"))
  })
})

test_that("comet", {
  set.seed(42)
  data("mtcars")
  expect_no_error({
    comet(mpg ~ cyl | disp, data = mtcars)
    comet(cbind(mpg, cyl) ~ vs + mpg | disp, data = mtcars)
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
