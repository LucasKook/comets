test_that("max GCM is standardized and invariant to rescaling", {
  set.seed(1)
  x <- rnorm(100)
  y <- rnorm(100)
  set.seed(2)
  base <- rgcm(y, x, type = "max", B = 199)
  set.seed(2)
  scaled <- rgcm(y, 100 * x, type = "max", B = 199)
  expect_equal(base$statistic, scaled$statistic)
  expect_equal(base$p.value, scaled$p.value)
  expect_equal(unname(base$statistic),
    abs(sqrt(length(x)) * mean(x * y) / sqrt(mean((x * y - mean(x * y))^2))))
  set.seed(2)
  greater <- rgcm(y, x, type = "max", alternative = "greater", B = 199)
  set.seed(2)
  less <- rgcm(y, -x, type = "max", alternative = "less", B = 199)
  expect_equal(unname(greater$statistic), unname(less$statistic))
  expect_equal(greater$p.value, less$p.value)
  expect_error(rgcm(rep(1, 10), rep(1, 10), type = "max"), "variance")
})

test_that("wGCM without coin uses the one-sided max statistic", {
  set.seed(43)
  n <- 100
  z <- matrix(rnorm(n), ncol = 1)
  x <- rnorm(n); y <- x * (z[, 1] > 0) + rnorm(n)
  set.seed(44)
  fit <- wgcm(y, x, z, reg_YonZ = "lrm", reg_XonZ = "lrm", reg_wfun = "lrm",
              B = 99, coin = FALSE)
  set.seed(44)
  sample.int(n, ceiling(n / 2))
  reference <- comets:::.gcm(fit$rY, as.matrix(fit$rX * fit$W),
                            alternative = "greater", type = "max", B = 99)
  expect_equal(unname(fit$statistic), reference$stat)
  expect_equal(fit$p.value, reference$pval)
})

test_that("quantile forest grouped residuals match the pairwise definition", {
  set.seed(91)
  x <- matrix(rnorm(80), 40, 2, dimnames = list(NULL, c("x1", "x2")))
  y <- sample(1:5, 40, replace = TRUE)
  model <- comets:::qrf(y, x, num.trees = 15, num.threads = 1, min.node.size = 2)
  nodes <- predict(model, data = x, type = "terminalNodes")$predictions
  k <- matrix(0, 40, 40)
  for (tree in seq_len(ncol(nodes))) k <- k + outer(nodes[, tree], nodes[, tree], "==")
  diag(k) <- 0
  k <- sweep(k, 1, pmax(colSums(k), .Machine$double.eps), "/")
  expected <- 2 * rowSums(k * outer(y, y, ">=")) - 1
  expect_equal(unname(residuals(model, data = x)), unname(expected), ignore_attr = TRUE)
  expect_error(residuals(model, data = x[-1, , drop = FALSE]), "training observations")
  expect_error(residuals(model, data = x[40:1, ]), "original order")
})

test_that("tuned forests retain the best OOB fit", {
  set.seed(10)
  x <- matrix(rnorm(160), 80, 2); y <- x[, 1] + rnorm(80)
  colnames(x) <- c("x1", "x2")
  set.seed(11)
  candidates <- lapply(1:2, function(depth) ranger::ranger(y = y, x = x, max.depth = depth, mtry = 1, num.trees = 20, num.threads = 1))
  set.seed(11)
  selected <- comets:::tuned_rf(y, x, max.depths = 1:2, mtrys = list(1), num.trees = 20, num.threads = 1)
  expect_equal(selected$prediction.error, min(vapply(candidates, `[[`, numeric(1), "prediction.error")))
})

test_that("tuned forests evaluate each distinct mtry once per depth", {
  set.seed(12)
  x <- matrix(rnorm(160), 80, 2, dimnames = list(NULL, c("x1", "x2")))
  y <- rnorm(80)
  options <- list(max.depths = 1:2, num.trees = 10, num.threads = 1, verbose = TRUE)
  set.seed(13)
  output <- capture.output(fit <- do.call(comets:::tuned_rf, c(list(y = y, x = x), options)))
  expect_length(output, 4L)
  set.seed(13)
  capture.output(reference <- do.call(comets:::tuned_rf, c(list(y = y, x = x, mtrys = list(1, 2)), options)))
  expect_equal(fit$prediction.error, reference$prediction.error)
  expect_equal(predict(fit, data = x), predict(reference, data = x))
  one_column <- capture.output(single_fit <- do.call(comets:::tuned_rf,
    c(list(y = y, x = x[, 1, drop = FALSE]), options)))
  expect_length(one_column, 2L)
})

test_that("tuned boosting refits at the best cross-validation iteration", {
  skip_if_not_installed("xgboost")
  set.seed(15)
  x <- matrix(rnorm(180), 60, 3); y <- rnorm(60)
  folds <- split(seq_len(60), rep(1:3, length.out = 60))
  params <- list(eta = 0.1, max_depth = 1, nthread = 1)
  set.seed(17)
  cv <- xgboost::xgb.cv(params = params, data = xgboost::xgb.DMatrix(x, label = y),
                        nrounds = 5, folds = folds, metrics = list("rmse"), verbose = 0)
  set.seed(17)
  model <- comets:::tuned_xgb(y, x, folds = folds, etas = 0.1, max_depths = 1,
                              nrounds = c(2, 5), nthread = 1)
  expect_equal(xgboost::xgb.get.num.boosted.rounds(model), which.min(cv$evaluation_log$test_rmse_mean))
})

test_that("kGCM batches the exact bootstrap and handles bandwidths", {
  set.seed(51)
  n <- 50
  y <- rnorm(n); x <- cbind(rnorm(n), rnorm(n)); z <- matrix(rnorm(n), ncol = 1)
  set.seed(7)
  fit <- kgcm(y, x, z, reg_YonZ = "lrm", reg_XonZ = "lrm", B = 23, batch_size = 4)
  rp <- c(fit$rY) * fit$rX
  set.seed(7)
  samples <- replicate(23, {
    r <- rp * sample(c(-1, 1), n, TRUE)
    sum(r * (fit$K %*% r)) / n
  })
  expect_equal(fit$p.value, (1 + sum(samples >= fit$statistic)) / 24)
  set.seed(7)
  compact <- kgcm(y, x, z, reg_YonZ = "lrm", reg_XonZ = "lrm", B = 23,
                  batch_size = 1, return_kernel = FALSE)
  expect_equal(compact$p.value, fit$p.value)
  expect_null(compact$K)
  expect_error(plot(compact), "not retained")
  expect_error(kgcm(y, x, z, bandwidth = 0, reg_YonZ = "lrm", reg_XonZ = "lrm"), "positive")
  expect_error(kgcm(y, x, z, B = 0), "positive integer")
  constant <- kgcm(y, x, matrix(1, n, 1), reg_YonZ = "lrm", reg_XonZ = "lrm", B = 9)
  expect_true(is.finite(constant$p.value))
  expect_true(all(constant$K == 1))
})

test_that("PCM caches variance predictions during calibration", {
  set.seed(32)
  n <- 100
  z <- matrix(rnorm(n), ncol = 1)
  x <- matrix(rnorm(n), ncol = 1)
  y <- x[, 1] + z[, 1] + rnorm(n)
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fit_variance <- function(y, x) {
    structure(list(value = mean(y), calls = calls), class = "counted_pcm_variance")
  }
  predict.counted_pcm_variance <- function(object, data, ...) {
    object$calls$n <- object$calls$n + 1L
    rep(object$value, NROW(data))
  }
  registerS3method("predict", "counted_pcm_variance", predict.counted_pcm_variance)
  on.exit(rm("predict.counted_pcm_variance", envir = get(".__S3MethodsTable__.", asNamespace("stats"))), add = TRUE)
  fit <- pcm(y, x, z, indices = 1:50, reg_YonXZ = "lrm", reg_YhatonZ = "lrm",
             reg_VonXZ = fit_variance, reg_YonZ = "lrm", reg_RonZ = "lrm")
  expect_true(is.finite(fit$p.value))
  expect_equal(calls$n, 2L)
})

test_that("entry points reject invalid observations and splits", {
  set.seed(4)
  y <- rnorm(20); x <- rnorm(20); z <- rnorm(20)
  expect_error(rgcm(y, x[-1]), "same number")
  expect_error(rgcm(y, replace(x, 1, Inf)), "non-finite")
  expect_error(rgcm(y, x, type = "max", B = 0), "positive integer")
  expect_equal(rgcm(y, data.frame(x = x))$p.value, rgcm(y, x)$p.value)
  expect_error(pcm(y, x[-1], z), "same number")
  expect_error(wgcm(y, x, replace(z, 1, NA)), "missing")
  expect_error(pcm(y, x, z, frac = 1), "frac")
  expect_error(wgcm(y, x, z, frac = 0.99), "at least two")
  expect_error(pcm(y, x, z, indices = c(1, 1, 2)), "unique")
  expect_error(pcm(y, x, z, indices = 1:21), "valid")
  expect_error(pcm(y, x, z, rep = 2, indices = 1:10), "list")
  expect_error(wgcm(y, x, z, B = 0), "positive integer")
})

test_that("quadratic GCM uses the rank of the covariance", {
  set.seed(3)
  y <- rnorm(120)
  x <- rnorm(120)
  scalar <- rgcm(y, x)
  expect_warning(duplicated <- rgcm(y, cbind(x, x)), "rank deficient")
  expect_equal(unname(duplicated$statistic), unname(scalar$statistic))
  expect_equal(unname(duplicated$parameter), 1)
  expect_equal(duplicated$p.value, scalar$p.value)
  x2 <- cbind(x, rnorm(120))
  rr <- y * x2
  covariance <- crossprod(scale(rr, scale = FALSE)) / length(y)
  expected <- length(y) * crossprod(colMeans(rr), solve(covariance, colMeans(rr)))
  full <- rgcm(y, x2)
  expect_equal(unname(full$statistic), as.numeric(expected))
  expect_equal(unname(full$parameter), 2)
  expect_error(rgcm(rep(1, 10), rep(1, 10)), "variance")
  expect_error(rgcm(rep(1, 10), matrix(1, 10, 2)), "zero rank")
})
