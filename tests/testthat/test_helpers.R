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
    gcm2 <- gcm(Y, X, Z,
      reg_XonZ = "lasso", reg_YonZ = "rf",
      args_YonZ = list(mtry = 2)
    )
    gcm3 <- gcm(Y, X, Z, reg_XonZ = "ridge", reg_YonZ = "ridge")
    gcm4 <- gcm(Y, X, Z, reg_XonZ = "qrf", reg_YonZ = "qrf")
    gcm5 <- gcm(Y, X, Z, reg_XonZ = "postlasso", reg_YonZ = "postlasso")
    gcm6 <- gcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm")
    if (require("xgboost")) {
      gcm7 <- gcm(Y, X, Z, reg_XonZ = "xgb", reg_YonZ = "xgb")
    }
    tmp <- plot(gcm6, plot = FALSE)
    gcm8 <- gcm(Y, X, Z, reg_XonZ = "unilasso", reg_YonZ = "unilasso")
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
  expect_no_error(comet(Surv(time, cens) ~ horTh | age,
    data = GBSG2,
    reg_YonZ = "cox"
  ))
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
      distribution = approximate(99)
    ))
    gcm(Y, X, Z, type = "quadratic", coin = TRUE, cointrol = list(
      distribution = "asymptotic"
    ))
    gcm(Y, X, Z, type = "max", coin = TRUE, cointrol = list(
      distribution = "asymptotic"
    ), alternative = "less")
    gcm(Y, X, Z, type = "max", coin = TRUE, cointrol = list(
      distribution = "asymptotic"
    ), alternative = "greater")
    gcm(Y[, 1], X[, 1], Z, type = "scalar", coin = TRUE, cointrol = list(
      distribution = "asymptotic"
    ))
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

test_that("comet wrapper works", {
  set.seed(42)
  data("mtcars")
  expect_no_error({
    comet(mpg ~ cyl | disp, data = mtcars)
    comet(mpg ~ cyl | 1, data = mtcars)
    comet(mpg ~ cyl | 1, data = mtcars, test = "pcm")
    comet(mpg ~ cyl | 1, data = mtcars, test = "wgcm")
    comet(cbind(mpg, cyl) ~ vs + mpg | disp, data = mtcars)
    comet(factor(vs) ~ cyl | disp, data = mtcars)
    comet(mpg ~ cyl | disp, data = mtcars, test = "pcm")
    comet(factor(vs) ~ cyl | disp, data = mtcars, test = "pcm")
    comet(mpg ~ cyl | disp, data = mtcars, test = "wgcm")
    comet(factor(vs) ~ cyl | disp, data = mtcars, test = "wgcm")
    comet(factor(mpg) ~ cyl | disp, data = mtcars)
    comet(ordered(mpg) ~ cyl | disp, data = mtcars)
    comet(cbind(mpg, vs) ~ cyl | disp, data = mtcars)
  })
  expect_error({
    comet(factor(mpg) ~ cyl | disp, data = mtcars, test = "pcm")
  })
  expect_error({
    comet(factor(mpg) ~ cyl | disp, data = mtcars, test = "wgcm")
  })
})

test_that("residual gcm works", {
  set.seed(42)
  expect_no_error({
    rY <- rnorm(100)
    rX <- rnorm(100)
    rX2 <- rnorm(100)
    rgcm(rY, rX, type = "quadratic")
    rgcm(rY, rX, type = "max")
    rgcm(rY, rX, type = "scalar")
    rgcm(rY, cbind(rX, rX2), type = "quadratic")
    rgcm(cbind(rY, rX2), rX, type = "max")
  })
})

test_that("glm regressions work", {
  m <- glrm(rnorm(10), rnorm(10))
  expect_length(residuals(m, response = rnorm(5), data = rnorm(5)), 5)
  expect_no_error({
    comet(Sepal.Width ~ Sepal.Length | Species,
      data = iris, reg_YonZ = "glrm",
      args_YonZ = list(family = "quasipoisson"),
      reg_XonZ = "glrm", args_XonZ = list(family = "Gamma")
    )
  })
  expect_no_error({
    comet(Sepal.Width ~ Sepal.Length | Species,
      test = "pcm",
      data = iris, reg_YonXZ = "glrm",
      args_YonXZ = list(family = "quasipoisson"),
      reg_YonZ = "glrm", args_YonZ = list(family = "Gamma")
    )
  })
})

test_that("rf works with different response types", {
  set.seed(12)
  tn <- 3e2
  dat <- data.frame(
    bin = factor(sample(0:1, tn, TRUE)),
    ord = ordered(sample(7:10, tn, TRUE)),
    mcc = factor(sample(11:15, tn, TRUE)),
    num = rnorm(tn),
    x = rnorm(tn)
  )
  lapply(colnames(dat)[-ncol(dat)], \(resp) {
    rf <- rf(y = dat[[resp]], x = dat[, "x", drop = FALSE])
    rr <- residuals.rf(rf,
      data = dat[, "x", drop = FALSE],
      response = dat[[resp]]
    )
    expect_true(all(colMeans(as.matrix(rr)) <= 0.05))
  })
  expect_no_error({
    comet(mcc ~ bin | x, data = dat)
  })
})

test_that("GCM with multivariate regressions works", {
  set.seed(12)
  tn <- 3e2
  X <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(X) <- c("X1", "X2")
  Z <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- cbind(Y1 = rnorm(tn), Y2 = rnorm(tn), Y3 = rnorm(tn))
  expect_no_error(gcm1 <- gcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm", multivariate = "both"))
})

test_that("GCM with (tuned) xgboost and rangers works", {
  if (require("xgboost") & require("lightgbm")) {
    set.seed(12)
    tn <- 1e2
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    expect_no_error(gcm(Y, X, Z, reg_XonZ = "tuned_xgb", reg_YonZ = "tuned_xgb"))
    expect_no_error(gcm(Y, X, Z, reg_XonZ = "xgb", reg_YonZ = "xgb"))
    expect_no_error(gcm(Y, X, Z, reg_XonZ = "tuned_rf", reg_YonZ = "tuned_rf"))
    expect_no_error(gcm(Y, X, Z, reg_XonZ = "lgbm", reg_YonZ = "lgbm"))
  }
})
