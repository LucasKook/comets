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
    gcm3 <- gcm(Y, as.data.frame(X), as.data.frame(Z))
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

test_that("kGCM works", {
  expect_no_error({
    set.seed(12)
    tn <- 1e2
    gamma <- 0.5
    # Same setup as simulation in original kgcm paper
    Z <- rnorm(tn)
    X <- Z + (U1 <- rnorm(tn)) * sin(5 * Z)
    Y <- Z^2 + gamma * U1 + (1 - gamma) * rnorm(tn)
    kgcm1 <- kgcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm")
    kgcm2 <- kgcm(Y, X, Z, bandwidth = 0.1)
    kgcm3 <- kgcm(Y, X, Z, bandwidth = 0.1, B = 999)
    tmp <- plot(kgcm3, plot = FALSE)
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

test_that("original for tests", {
  set.seed(12)
  tn <- 3e2
  X <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(X) <- c("X1", "X2")
  Z <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- cbind(rowSums(X) + rnorm(tn), rowSums(X) + rnorm(tn))
  expect_no_error({
    gcm(Y, X, Z, type = "max", coin = FALSE)
    gcm(Y, X, Z, type = "quadratic", coin = FALSE)
    gcm(Y[, 1], X, Z, type = "quadratic", coin = FALSE)
    gcm(Y[, 1], X, Z, type = "max", coin = FALSE)
    gcm(Y[, 1], X[, 1], Z, type = "quadratic", coin = FALSE)
    gcm(Y[, 1], X[, 1], Z, type = "max", coin = FALSE)
  })
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
    comet(factor(cyl) ~ mpg | disp, data = mtcars)
    comet(ordered(cyl) ~ mpg | disp, data = mtcars)
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
    rgcm(rY, rX, type = "quadratic", coin = FALSE)
    rgcm(rY, rX, type = "max", coin = FALSE)
    rgcm(rY, rX, type = "scalar", coin = FALSE)
    rgcm(rY, cbind(rX, rX2), type = "quadratic", coin = FALSE)
    rgcm(cbind(rY, rX2), rX, type = "max", coin = FALSE)
    rgcm(rY, rX, type = "quadratic", coin = TRUE)
    rgcm(rY, rX, type = "max", coin = TRUE)
    rgcm(rY, rX, type = "scalar", coin = TRUE)
    rgcm(rY, cbind(rX, rX2), type = "quadratic", coin = TRUE)
    rgcm(cbind(rY, rX2), rX, type = "max", coin = TRUE)
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

test_that("rf residuals for ordered responses work", {
  set.seed(12)
  tn <- 2e2
  z <- rnorm(tn)
  x <- matrix(z, dimnames = list(NULL, "z"))
  y <- cut(z + rlogis(tn), c(-Inf, -1, 0, 1.5, Inf), ordered_result = TRUE)
  ### Closed form score residual agrees with MASS::polr
  if (requireNamespace("MASS") && requireNamespace("sandwich")) {
    m <- MASS::polr(y ~ z)
    sc <- sandwich::estfun(m)
    r_polr <- -unname(rowSums(sc[, names(m$zeta)]))
    expect_equal(.ordinal_residuals(y, stats::fitted(m), "score"), r_polr)
  }
  ### Dimensions and mean zero
  m <- rf(y, x)
  expect_equal(dim(residuals(m, response = y, data = x)), c(tn, 3L))
  for (type in c("mean", "score")) {
    m <- rf(y, x, ordinal_residuals = type)
    r <- residuals(m, response = y, data = x)
    expect_length(r, tn)
    expect_lt(abs(mean(r)), 0.1)
    m <- tuned_rf(y, x, ordinal_residuals = type, max.depths = 1:2)
    expect_length(residuals(m, response = y, data = x), tn)
  }
  ### Score residuals are bounded, and mean residuals use integer scores
  probs <- matrix(c(0.2, 0.3, 0.5), nrow = 3, ncol = 3, byrow = TRUE)
  yy <- ordered(1:3)
  expect_equal(.ordinal_residuals(yy, probs, "score"), c(-0.8, -0.3, 0.5))
  expect_equal(.ordinal_residuals(yy, probs, "mean"), 1:3 - 2.3)
  ### Binary ordered responses are unaffected
  yb <- ordered(as.numeric(y) > 2)
  expect_equal(
    residuals(rf(yb, x, ordinal_residuals = "score", seed = 1),
      response = yb, data = x
    ),
    residuals(rf(yb, x, seed = 1), response = yb, data = x)
  )
  ### Usage through comet()
  dat <- data.frame(y = y, z = z, x = rnorm(tn))
  for (type in c("probability", "mean", "score")) {
    expect_no_error(tst <- comet(y ~ x | z,
      data = dat,
      args_YonZ = list(ordinal_residuals = type)
    ))
    expect_equal(NCOL(tst$rY), if (type == "probability") 3L else 1L)
  }
})

test_that("rf drops unused response levels with a warning", {
  set.seed(12)
  tn <- 1e2
  lv <- c("none", "mild", "moderate", "severe")
  d <- data.frame(
    Y = ordered(sample(lv[-2], tn, TRUE), levels = lv),
    X = rnorm(tn), Z = rnorm(tn)
  )
  x <- as.matrix(d["Z"])
  for (type in c("probability", "mean", "score")) {
    expect_warning(
      m <- rf(d$Y, x, ordinal_residuals = type),
      "Dropping unused level(s) of the response: mild",
      fixed = TRUE
    )
    r <- residuals(m, response = d$Y, data = x)
    expect_equal(NCOL(r), if (type == "probability") 2L else 1L)
  }
  expect_warning(m <- tuned_rf(d$Y, x, max.depths = 1:2), "mild")
  expect_equal(NCOL(residuals(m, response = d$Y, data = x)), 2L)
  expect_warning(tst <- comet(Y ~ X | Z, data = d), "mild")
  expect_equal(NCOL(tst$rY), 2L)
  ### Two remaining levels are treated as binary
  d2 <- subset(d, Y != "moderate")
  expect_warning(m <- rf(d2$Y, as.matrix(d2["Z"])), "mild, moderate")
  expect_length(residuals(m, response = d2$Y, data = as.matrix(d2["Z"])), NROW(d2))
  ### Evaluation data missing a level seen during training
  y <- ordered(sample(lv, tn, TRUE), levels = lv)
  m <- rf(y, x, ordinal_residuals = "score")
  idx <- which(y != "mild")
  expect_length(residuals(m, response = y[idx], data = x[idx, , drop = FALSE]), length(idx))
  m <- rf(y, x)
  expect_equal(dim(residuals(m, response = y[idx], data = x[idx, , drop = FALSE])), c(length(idx), 3L))
  ### Evaluation data containing a level not seen during training
  expect_warning(m <- rf(d$Y, x))
  expect_error(residuals(m, response = y, data = x), "not seen during training")
})

test_that("printed hypothesis and method depend on the residuals", {
  set.seed(12)
  tn <- 1e2
  d <- data.frame(Z = rnorm(tn), X = rnorm(tn), Yn = rnorm(tn))
  d$Y <- cut(d$Z + rlogis(tn), c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
  hyp <- \(tst) names(tst$null.value)
  tst <- comet(Yn ~ X | Z, data = d)
  expect_equal(hyp(tst), "E[cov(Y, X | Z)]")
  expect_equal(tst$method, "Generalised covariance measure test")
  expect_equal(hyp(comet(Y ~ X | Z, data = d)), "E[cov(1(Y = k), X | Z)]")
  expect_equal(hyp(comet(Y ~ X | Z,
    data = d,
    args_YonZ = list(ordinal_residuals = "mean")
  )), "E[cov(Y, X | Z)]")
  for (coin in c(FALSE, TRUE)) {
    tst <- comet(Y ~ X | Z,
      data = d, coin = coin,
      args_YonZ = list(ordinal_residuals = "score")
    )
    expect_equal(hyp(tst), "E[cov(R, X | Z)]")
    expect_equal(tst$method, "TRAM-GCM test")
  }
  expect_equal(
    hyp(comet(Yn ~ X | Z, data = d, reg_YonZ = "qrf")),
    "E[cov(F(Y | Z), X | Z)]"
  )
  expect_equal(
    hyp(comet(Yn ~ X | Z, data = d, test = "wgcm", reg_YonZ = "qrf")),
    "E[w(Z) cov(F(Y | Z), X | Z)]"
  )
  expect_equal(hyp(kgcm(d$Yn, d$X, d$Z)), "E[w(Z) cov(Y, X | Z)]")
  library("survival")
  data("GBSG2", package = "TH.data")
  for (reg in c("cox", "survforest")) {
    tst <- comet(Surv(time, cens) ~ horTh | age, data = GBSG2, reg_YonZ = reg)
    expect_equal(hyp(tst), "E[cov(R, X | Z)]")
    expect_equal(tst$method, "TRAM-GCM test")
  }
  ### Pre-computed residuals
  expect_equal(hyp(rgcm(rnorm(tn), rnorm(tn))), "E[cov(Y, X | Z)]")
  expect_equal(
    hyp(rgcm(rnorm(tn), rnorm(tn), hypothesis = "E[cov(R, X | Z)]")),
    "E[cov(R, X | Z)]"
  )
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
    expect_no_error(gcm(Y, X, Z,
      reg_XonZ = "rf", reg_YonZ = "tuned_xgb",
      args_YonZ = list(nfold = 2)
    ))
    expect_no_error(gcm(Y, X, Z,
      reg_XonZ = "rf", reg_YonZ = "tuned_xgb",
      args_YonZ = list(folds = list(1:50, 51:100))
    ))
  }
})

test_that("PCM with rep larger 1 and indices", {
  expect_no_error({
    set.seed(12)
    tn <- 3e2
    set.seed(12)
    X <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2 * tn), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(tn)
    pcm(Y, X, Z, reg_YonXZ = "lrm", reg_YonZ = "lrm", rep = 2)
    rm <- "lrm"
    pcm(Y, X, Z, reg_YonXZ = rm, reg_YonZ = rm, rep = 2)
    indices <- list(1:150, 151:300)
    pcm(Y, X, Z, reg_YonXZ = "lrm", reg_YonZ = "lrm", rep = 2, indices = indices)
    tmp <- \(rm) {
      pcm(Y, X, Z, reg_YonXZ = rm, reg_YonZ = rm, rep = 2)
    }
    tmp("lrm")
  })
})

test_that("kGCM works with multivariate X and matrix/data.frame inputs", {
  set.seed(12)
  tn <- 1e2
  X <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(X) <- c("X1", "X2")
  Z <- matrix(rnorm(2 * tn), ncol = 2)
  colnames(Z) <- c("Z1", "Z2")
  Y <- rnorm(tn)
  kgcm1 <- kgcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm", B = 99)
  expect_length(kgcm1$statistic, 1)
  expect_true(kgcm1$p.value > 0 && kgcm1$p.value <= 1)
  expect_no_error(kgcm(cbind(Y), X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm", B = 99))
  expect_no_error(kgcm(Y, X, as.data.frame(Z), reg_XonZ = "lrm", reg_YonZ = "lrm", B = 99))
  ### Univariate X gives the same statistic as before
  kgcm2 <- kgcm(Y, X[, 1], Z, reg_XonZ = "lrm", reg_YonZ = "lrm", B = 99)
  RP <- c(kgcm2$rY) * c(kgcm2$rX)
  expect_equal(unname(kgcm2$statistic), c(t(RP) %*% kgcm2$K %*% RP) / tn)
})

test_that("scalar GCM returns signed statistic", {
  set.seed(42)
  rY <- rnorm(100)
  rX <- rnorm(100) - rY
  tst <- rgcm(rY, rX, type = "scalar", coin = FALSE)
  expect_true(tst$statistic < 0)
  expect_null(tst$parameter)
  expect_equal(tst$p.value, 2 * pnorm(-abs(unname(tst$statistic))))
  expect_equal(
    unname(rgcm(rY, rX, type = "quadratic", coin = FALSE)$statistic),
    unname(tst$statistic)^2
  )
  expect_error(rgcm(rY, cbind(rX, rnorm(100)), type = "scalar", coin = FALSE))
  expect_error(rgcm(rY, rX, alternative = "foo"))
})

test_that("data.frame inputs and missing values are handled", {
  set.seed(12)
  tn <- 1e2
  X <- data.frame(X1 = rnorm(tn), X2 = rnorm(tn))
  Z <- data.frame(Z1 = rnorm(tn), Z2 = rnorm(tn))
  Y <- data.frame(Y1 = rnorm(tn), Y2 = rnorm(tn))
  expect_no_error(gcm(Y, X, Z, reg_XonZ = "lrm", reg_YonZ = "lrm"))
  expect_no_error(gcm(Y[, 1], X, Z, reg_XonZ = "postlasso", reg_YonZ = "glrm"))
  expect_no_error(pcm(Y[, 1], X, Z, reg_YonXZ = "lrm", reg_YonZ = "lrm"))
  dat <- data.frame(y = Y[, 1], x = X[, 1], z = Z[, 1])
  dat$y[1:5] <- NA
  dat$z[6:10] <- NA
  expect_no_error(tst <- comet(y ~ x | z, data = dat, reg_XonZ = "lrm", reg_YonZ = "lrm"))
  expect_length(tst$rY, tn - 10)
})

test_that("qrf residuals are conditional on Z", {
  set.seed(12)
  tn <- 2e2
  z <- rnorm(tn)
  y <- 3 * z + rnorm(tn)
  m <- qrf(y, cbind(z))
  r <- residuals(m, response = y, data = cbind(z))
  expect_true(abs(cor(r, z)) < 0.2)
  expect_true(all(abs(r) <= 1 + 1e-8))
})

test_that("wGCM passes cointrol to coin", {
  set.seed(12)
  tn <- 1e2
  X <- rnorm(tn)
  Z <- rnorm(tn)
  Y <- rnorm(tn)
  tst <- wgcm(Y, X, Z,
    reg_XonZ = "lrm", reg_YonZ = "lrm", reg_wfun = "lrm",
    cointrol = list(distribution = coin::approximate(9))
  )
  expect_true(any(abs(tst$p.value - (0:9) / 9) < 1e-8))
})
