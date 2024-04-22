
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
    gcm2 <- gcm(Y, X, Z, reg_XonZ = "lasso", reg_YonZ = "rf")
    gcm3 <- gcm(Y, X, Z, reg_XonZ = "ridge", reg_YonZ = "ridge")
    gcm4 <- gcm(Y, X, Z, reg_XonZ = "qrf", reg_YonZ = "qrf")
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
  debugonce(comet)
  expect_no_error(comet(Surv(time, cens) ~ horTh | age, data = GBSG2,
                        reg_YonZ = "cox"))
})
