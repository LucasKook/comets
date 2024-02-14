
test_that(".ranger works", {
  set.seed(12)
  dat <- data.frame(bin = factor(sample(0:1, 1e3, TRUE)),
                    ord = ordered(sample(7:10, 1e3, TRUE)),
                    mcc = factor(sample(11:15, 1e3, TRUE)),
                    num = rnorm(1e3),
                    x = rnorm(1e3))
  lapply(colnames(dat)[-ncol(dat)], \(resp) {
    fm <- reformulate("x", resp)
    rf <- .ranger(fm, data = dat)
    rr <- residuals.ranger(rf)
    expect_lt(abs(mean(rr)), 0.005)
  })
})

test_that("gcm and pcm work", {
  expect_no_error({
    set.seed(12)
    X <- matrix(rnorm(2e3), ncol = 2)
    colnames(X) <- c("X1", "X2")
    Z <- matrix(rnorm(2e3), ncol = 2)
    colnames(Z) <- c("Z1", "Z2")
    Y <- rnorm(1e3)
    gcm1 <- gcm(Y, X, Z)
    pcm1 <- pcm(Y, X, Z)
  })
})
