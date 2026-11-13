experimental_equiv <- function(...) suppressWarnings(plm_equiv_test(...))

test_that("experimental equivalence scales warn", {
  set.seed(417)
  x <- rnorm(50); y <- rnorm(50); z <- matrix(rnorm(50), ncol = 1)
  args <- list(Y = y, X = x, Z = z, from = -1, to = 1,
               reg_YonZ = "lrm", reg_XonZ = "lrm")
  expect_warning(do.call(plm_equiv_test, args), "'plm'.*experimental")
  expect_warning(do.call(plm_equiv_test, c(args, list(scale = "cor"))), "'cor'.*experimental")
  expect_no_warning(do.call(plm_equiv_test, c(args, list(scale = "cov"))))
})

test_that("ratio equivalence tests use their influence functions", {
  set.seed(418)
  n <- 120
  z <- matrix(rnorm(n), ncol = 1)
  x <- 0.4 * z[, 1] + rnorm(n)
  y <- 0.8 * x + z[, 1] + rnorm(n)
  for (scale in c("plm", "cor")) {
    fit <- experimental_equiv(y, x, z, from = 0.1, to = 1.1, scale = scale,
                          reg_YonZ = "lrm", reg_XonZ = "lrm", type = "scalar")
    u <- c(fit$rX); v <- c(fit$rY)
    xx <- mean(u^2); yy <- mean(v^2); xy <- mean(u * v)
    estimate <- if (scale == "plm") xy / xx else xy / sqrt(xx * yy)
    influence <- if (scale == "plm") {
      u * (v - estimate * u) / xx
    } else {
      u * v / sqrt(xx * yy) - estimate / 2 * (u^2 / xx + v^2 / yy)
    }
    se <- sqrt(var(influence) / n)
    expected <- pchisq(((estimate - 0.6) / se)^2, df = 1,
                       ncp = (0.5 / se)^2)
    expect_equal(unname(fit$estimate), estimate)
    expect_equal(fit$std.error, se)
    expect_equal(fit$p.value, expected)
  }
})

test_that("equivalence inference is invariant to positive changes of units", {
  set.seed(419)
  n <- 120
  z <- matrix(rnorm(n), ncol = 1)
  x <- rnorm(n); y <- 0.5 * x + rnorm(n)
  opts <- list(reg_YonZ = "lrm", reg_XonZ = "lrm")
  cor1 <- do.call(experimental_equiv, c(list(y, x, z, from = -0.3, to = 0.8, scale = "cor"), opts))
  cor2 <- do.call(experimental_equiv, c(list(10 * y, 3 * x, z, from = -0.3, to = 0.8, scale = "cor"), opts))
  expect_equal(cor1$p.value, cor2$p.value)
  plm1 <- do.call(experimental_equiv, c(list(y, x, z, from = -0.3, to = 0.8, scale = "plm"), opts))
  plm2 <- do.call(experimental_equiv, c(list(10 * y, 3 * x, z, from = -1, to = 8 / 3, scale = "plm"), opts))
  expect_equal(plm1$p.value, plm2$p.value)
})

test_that("covariance equivalence retains its original statistic", {
  set.seed(421)
  n <- 120
  z <- matrix(rnorm(n), ncol = 1)
  x <- rnorm(n); y <- 0.2 * x + rnorm(n)
  fit <- plm_equiv_test(y, x, z, from = -0.2, to = 0.5, scale = "cov",
                        reg_YonZ = "lrm", reg_XonZ = "lrm")
  products <- c(fit$rX) * c(fit$rY)
  variance <- var(products)
  expected <- pchisq(n * (mean(products) - 0.15)^2 / variance,
                     df = 1, ncp = n * 0.35^2 / variance)
  expect_equal(fit$p.value, expected)
})

test_that("ratio equivalence tests are calibrated at nonzero boundaries", {
  set.seed(420)
  nsim <- 1000L
  n <- 400L
  for (scale in c("plm", "cor")) {
    for (boundary in c(-0.5, 0.5)) {
      p <- replicate(nsim, {
        z <- matrix(rnorm(n), ncol = 1)
        x <- rnorm(n)
        y <- boundary * x + sqrt(1 - boundary^2) * rnorm(n) + z[, 1]
        experimental_equiv(y, x, z, from = -0.5, to = 0.5, scale = scale,
                       reg_YonZ = "lrm", reg_XonZ = "lrm")$p.value
      })
      # 0.02--0.08 is wider than four Monte Carlo standard errors at 0.05.
      expect_gt(mean(p < 0.05), 0.02)
      expect_lt(mean(p < 0.05), 0.08)
    }
  }
})
