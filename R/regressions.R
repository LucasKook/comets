
# Random forest -----------------------------------------------------------

#' Implemented regression methods
#' @rdname regressions
#' @param y Vector (or matrix) of response values.
#' @param x Design matrix of predictors.
#' @param ... Additional arguments passed to the underlying regression method.
#'    In case of \code{"rf"}, \code{"survforest"} and \code{"qrf"}, this is
#'    \code{\link[ranger]{ranger}}. In case of \code{"lasso"} and
#'    \code{"ridge"}, this is \code{\link[glmnet]{glmnet}}. In case of
#'    \code{"cox"}, this is \code{\link[survival]{coxph}}.
#' @details
#' The implemented choices are \code{"rf"} for random forests as implemented in
#' ranger, \code{"lasso"} for cross-validated Lasso regression, \code{"ridge"}
#' for cross-validated ridge regression, \code{"cox"} for the Cox proportional
#' hazards model as implemented in survival, \code{"qrf"} or \code{"survforest"}
#' for quantile and survival random forests, respectively.
#' New regression methods can be implemented and supplied as well and need the
#' following structure. The regression method \code{"custom_reg"} needs to take
#' arguments \code{y, x, ...}, fit the model using \code{y} and \code{x} as
#' matrices and return an object of a user-specified class, for instance,
#' '\code{custom}'. For the GCM test, implementing a \code{residuals.custom}
#' method is sufficient, which should take arguments
#' \code{object, response = NULL, data = NULL, ...}. For the PCM test, a
#' \code{predict.custom} method is necessary for out-of-sample prediction
#' and computation of residuals.
#'
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
#' @rdname regressions
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
#' @rdname regressions
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
#' @rdname regressions
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

#' @rdname regressions
ridge <- function(y, x, ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), alpha = 0, ...)
  class(obj) <- c("lasso", class(obj))
  obj
}

# Cox ---------------------------------------------------------------------

#' @rdname regressions
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
