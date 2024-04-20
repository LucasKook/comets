
# Random forest -----------------------------------------------------------

rf <- function(y, x, ...) {
  args <- list(...)
  if (length(unique(y)) == 2) {
    y <- as.factor(y)
    args$probability <- TRUE
  }
  rf <- do.call("ranger", c(list(y = y, x = x), args))
  class(rf) <- c("rf", class(rf))
  rf
}

#' @exportS3Method predict rf
predict.rf <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  preds <- predict(object, data = data)$predictions
  if (object$treetype == "Probability estimation")
    preds <- preds[, 2]
  preds
}

#' @exportS3Method residuals rf
residuals.rf <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.rf(object, data, ...)
  .compute_residuals(response, preds)
}

### Survival forest
survforest <- function(y, x, ...) {
  rf <- ranger::ranger(y = y, x = x, ...)
  class(rf) <- c("survforest", class(rf))
  rf
}

#' @exportS3Method residuals survforest
residuals.survforest <- function(object, response, data, ...) {
  times <- response[, 1]
  status <- response[, 2]
  pred <- stats::predict(object, data = data)
  idx <- sapply(times, \(x) which.min(abs(x - pred$unique.death.times))[1])
  preds <- pred$survival
  ipreds <- sapply(seq_len(nrow(preds)), \(smpl) {
    -log(preds[smpl, idx[smpl]])
  })
  .compute_residuals(status, ipreds)
}

### Quantile forest
qrf <- function(y, x, ...) {
  rf <- ranger::ranger(y = y, x = x, ...)
  rf$response <- y
  rf$x <- x
  class(rf) <- c("qrf", class(rf))
  rf
}

#' @exportS3Method residuals qrf
residuals.qrf <- \(object, data, ...) {
  class(object) <- class(object)[-1]
  tn <- stats::predict(object, data = data, type = "terminalNodes")$predictions
  K <- matrix(0, nrow = N <- nrow(tn), ncol = N)
  for (tree in seq_len(B <- object$num.trees)) {
    K <- K + sapply(seq_len(nrow(tn)), \(obs) {
      as.numeric(tn[obs, tree] == tn[, tree])
    }, simplify = "matrix")
  }
  K <- K / B
  diag(K) <- 0
  K <- K / pmax(colSums(K), .Machine$double.eps)
  pred <- \(y) mean(K %*% as.numeric(object$response <= y))
  ipreds <- sapply(object$response, pred)
  2 * ipreds - 1
}

# Lasso/ridge -------------------------------------------------------------

#' @importFrom glmnet cv.glmnet
lasso <- function(y, x, ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), ...)
  class(obj) <- c("lasso", class(obj))
  obj
}

#' @exportS3Method predict lasso
predict.lasso <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, newx = as.matrix(data), s = object$lambda.min)[, 1]
}

#' @exportS3Method residuals lasso
residuals.lasso <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.lasso(object, data)
  .compute_residuals(response, preds)
}

ridge <- function(y, x, ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), alpha = 0, ...)
  class(obj) <- c("lasso", class(obj))
  obj
}

# Cox ---------------------------------------------------------------------

cox <- function(y, x, ...) {
  obj <- survival::coxph(y ~ x, ...)
  class(obj) <- c("cox", class(obj))
  obj
}

#' @exportS3Method predict cox
predict.cox <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, newdata = as.data.frame(data), ...)
}

#' @exportS3Method residuals cox
residuals.cox <- function(object, response = NULL, data = NULL, ...) {
  preds <- -log(predict.cox(object, data, type = "survival"))
  .compute_residuals(response[, 2], preds)
}
