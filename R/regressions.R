#' Implemented regression methods
#' @rdname regressions
#' @param y Vector (or matrix) of response values.
#' @param x Design matrix of predictors.
#' @param ... Additional arguments passed to the underlying regression method.
#'    In case of \code{"rf"}, \code{"tuned_rf"}, \code{"survforest"} and
#'    \code{"qrf"}, this is \code{\link[ranger]{ranger}}. In case of
#'    \code{"lasso"} and \code{"ridge"}, this is \code{\link[glmnet]{glmnet}}.
#'    In case of \code{"cox"}, this is \code{\link[survival]{coxph}}. In case
#'    of \code{"xgb"} and \code{"tuned_xgb"} this is
#'    \code{\link[xgboost]{xgboost}}.
#' @details
#' The implemented choices are \code{"rf"} for random forests as implemented in
#' ranger, \code{"lasso"} for cross-validated Lasso regression (using the
#' one-standard error rule), \code{"ridge"}
#' for cross-validated ridge regression (using the one-standard error rule),
#' \code{"cox"} for the Cox proportional
#' hazards model as implemented in survival, \code{"qrf"} or \code{"survforest"}
#' for quantile and survival random forests, respectively. The option
#' \code{"postlasso"} option refers to a cross-validated LASSO (using the
#' one-standard error rule) and subsequent OLS regression. The \code{"lrm"}
#' option implements a standard linear regression model. The \code{"xgb"} and
#' \code{"tuned_xgb"} options require the \code{xgboost} package.
#'
#' The \code{"tuned_rf"} regression method tunes the \code{mtry} and
#' \code{max.depth} parameters in \code{\link[ranger]{ranger}} out-of-bag.
#' The \code{"tuned_xgb"} regression method uses k-fold cross-validation to
#' tune the \code{nrounds}, \code{mtry} and \code{max_depth} parameters in
#' \code{\link[xgboost]{xgb.cv}}.
#'
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
    y <- factor(y)
  }
  if (is.factor(y)) {
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
  if (object$treetype == "Probability estimation") {
    preds <- preds[, -1]
  }
  preds
}

#' @exportS3Method residuals rf
residuals.rf <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.rf(object, data, ...)
  if (length(unique(response)) == 2) {
    response <- factor(response)
  }
  .compute_residuals(response, preds)
}

### Internal residuals function dealing with factors
#' @importFrom stats model.matrix model.frame predict
.compute_residuals <- function(y, pred) {
  if (is.factor(y)) {
    y <- stats::model.matrix(~ 0 + y,
      contrasts.arg = list("y" = "contr.treatment")
    )[, -1]
  }
  y - pred
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

#' @rdname regressions
lrm <- function(y, x, ...) {
  obj <- stats::lm(y ~ x, ...)
  class(obj) <- c("lrm", class(obj))
  obj
}

#' @exportS3Method predict lrm
predict.lrm <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  cfx <- object$coefficients
  cfx[is.na(cfx)] <- 0
  c(cbind(1, as.matrix(data)) %*% cfx)
}

#' @exportS3Method residuals lrm
residuals.lrm <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.lrm(object, data)
  .compute_residuals(response, preds)
}

#' @rdname regressions
#' @importFrom stats glm
glrm <- function(y, x, ...) {
  dat <- list(y = y, x = x)
  obj <- stats::glm(y ~ x, data = dat, ...)
  class(obj) <- c("glrm", class(obj))
  obj
}

#' @exportS3Method predict glrm
predict.glrm <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, newdata = list(x = data), ...)
}

#' @exportS3Method residuals glrm
#' @importFrom stats residuals
residuals.glrm <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict(object, data = data, type = "response")
  .compute_residuals(response, preds)
}

#' @importFrom glmnet cv.glmnet
#' @param s Which lambda to use for prediction, defaults to
#'     \code{"lambda.min"}. See \code{\link[glmnet]{cv.glmnet}}
#' @rdname regressions
lasso <- function(y, x, s = "lambda.min", ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), ...)
  obj$s <- s
  class(obj) <- c("lasso", class(obj))
  obj
}

#' @exportS3Method predict lasso
predict.lasso <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, newx = as.matrix(data), s = object$s)[, 1]
}

#' @exportS3Method residuals lasso
residuals.lasso <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.lasso(object, data)
  .compute_residuals(response, preds)
}

#' @rdname regressions
ridge <- function(y, x, s = "lambda.min", ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), alpha = 0, ...)
  obj$s <- s
  class(obj) <- c("lasso", class(obj))
  obj
}

#' @rdname regressions
postlasso <- function(y, x, s = "lambda.min", ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), ...)
  nz <- which(stats::coef(obj, s = s)[-1] != 0)
  obj <- if (!identical(nz, integer(0))) {
    stats::lm(y ~ x[, nz])
  } else {
    stats::lm(y ~ 1)
  }
  obj$nz <- nz
  class(obj) <- c("postlasso", class(obj))
  obj
}

#' @exportS3Method predict postlasso
predict.postlasso <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  c(cbind(1, as.matrix(data)[, object[["nz"]]]) %*% object[["coefficients"]])
}

#' @exportS3Method residuals postlasso
residuals.postlasso <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict.postlasso(object, data)
  .compute_residuals(response, preds)
}

### unilasso
unilasso <- function(
    y, x, loo = TRUE, lower.limits = 0, standardize = FALSE,
    s = "lambda.min", ...) {
  ms <- apply(x, 2, \(tx) {
    stats::lm(y ~ tx)
  }, simplify = FALSE)
  xloo <- if (loo) {
    lapply(ms, \(m) {
      y - stats::residuals(m) / (1 - stats::hatvalues(m))
    })
  } else {
    lapply(ms, stats::predict)
  }
  xloo <- do.call("cbind", xloo)
  ml <- glmnet::cv.glmnet(
    x = xloo, y = y, standardize = standardize,
    lower.limits = lower.limits
  )
  cfs <- lapply(ms, stats::coef) |> do.call("rbind", args = _)
  cf0 <- cfs[, 1]
  cfx <- cfs[, 2]
  cfl <- as.double(cfll <- stats::coef(ml, s = s))
  coef <- c(cfl[1] + sum(cfl[-1] * cf0), cfl[-1] * cfx)
  names(coef) <- rownames(cfll)

  structure(list(
    # y = y,
    # x = x,
    # xloo = xloo,
    lm_fits = ms,
    lasso_fit = ml,
    cfs = cfs,
    cfl = cfl,
    coef = coef
  ), class = "unilasso")
}

#' @exportS3Method predict unilasso
predict.unilasso <- function(object, data, ...) {
  cbind(1, data) %*% object$coef
}

#' @exportS3Method residuals unilasso
residuals.unilasso <- function(object, response, data, ...) {
  response - stats::predict(object, data, ...)
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

# Tuned (mtry/max.depth) ranger ------------------------------------------

#' @rdname regressions
#' @param max.depths Values for \code{max.depth} to tune out-of-bag. See
#'     \code{\link[ranger]{ranger}}.
#' @param mtrys for \code{mtry} to tune out-of-bag. See
#'     \code{\link[ranger]{ranger}}.
tuned_rf <- function(y, x, max.depths = 1:5,
                     mtrys = list(1, \(p) ceiling(sqrt(p)), identity),
                     verbose = FALSE,
                     ...) {
  args <- list(...)
  if (length(unique(y)) == 2) {
    y <- factor(y)
  }
  if (is.factor(y)) {
    args$probability <- TRUE
  }

  ### Tune OOB max.depth
  tmp_args <- args
  rfs <- lapply(seq_along(max.depths), \(tmd) {
    lapply(seq_along(mtrys), \(tmt) {
      if (verbose) {
        cat(
          "Tuning step with max.depth", tmd, "out of", length(max.depths),
          "and mtry", tmt, "out of", length(mtrys), "\n"
        )
      }
      tmp_args$max.depth <- max.depths[tmd]
      tmp_args$mtry <- mtrys[[tmt]]
      do.call("ranger", c(list(y = y, x = x), tmp_args))
    })
  }) |> unlist(recursive = FALSE)

  woob <- which.min(sapply(rfs, `[[`, "prediction.error"))
  rf <- rfs[[woob]]
  class(rf) <- c("rf", class(rf))
  rf
}

# Boosting ---------------------------------------------------------------

#' @rdname regressions
#' @param nrounds See \code{\link[xgboost]{xgboost}}.
#' @param verbose See \code{\link[xgboost]{xgboost}}.
xgb <- function(y, x, nrounds = 2L, verbose = 0L, ...) {
  if (requireNamespace("xgboost")) {
    bst <- do.call("xgboost", c(list(
      data = x, label = y, nrounds = nrounds,
      verbose = verbose
    ), list(...)))
    class(bst) <- c("xgb", class(bst))
    return(bst)
  }
  stop("Package `xgboost` not available.")
}

#' @rdname regressions
#' @param etas Values for \code{eta} to cross-validate. See
#'     \code{\link[xgboost]{xgboost}}.
#' @param max_depths Values for \code{max_depth} to cross-validate. See
#'     \code{\link[xgboost]{xgboost}}.
#' @param metrics See \code{\link[xgboost]{xgboost}}.
#' @param nfold Number of folds for \code{nfold}-cross validation.
tuned_xgb <- function(y, x, etas = c(0.1, 0.5, 1), max_depths = 1:5,
                      nfold = 5, nrounds = c(2, 10, 50), verbose = 0,
                      metrics = list("rmse"), ...) {
  if (requireNamespace("xgboost")) {
    cvres <- lapply(etas, \(teta) {
      lapply(max_depths, \(tmd) {
        lapply(nrounds, \(tnr) {
          cv <- do.call("xgb.cv", c(list(
            data = x, label = y, nrounds = tnr,
            verbose = verbose, eta = teta, max_depth = tmd,
            metrics = metrics, nfold = nfold
          ), list(...)))
          err <- mean(cv$evaluation_log[[paste0("test_", metrics[[1]], "_mean")]])
          data.frame(
            nrounds = tnr, eta = teta, max_depth = tmd, error = err
          )
        }) |> do.call("rbind", args = _)
      }) |> do.call("rbind", args = _)
    }) |> do.call("rbind", args = _)
    best <- which.min(cvres$error)[1]
    bst <- xgb(y, x,
      nrounds = cvres$nrounds[best], verbose = verbose,
      max_depth = cvres$max_depth[best], eta = cvres$eta[best], ...
    )
    class(bst) <- c("xgb", class(bst))
    return(bst)
  }
  stop("Package `xgboost` not available.")
}

#' @exportS3Method predict xgb
predict.xgb <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, data, ...)
}

#' @exportS3Method residuals xgb
residuals.xgb <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict(object, data = data, ...)
  .compute_residuals(response, preds)
}

#' @rdname regressions
lgbm <- function(y, x, nrounds = 100L, verbose = -1L, ...) {
  if (requireNamespace("lightgbm")) {
    dtr <- lightgbm::lgb.Dataset(data = x, label = y)
    bst <- do.call("lgb.train", c(list(
      data = dtr, nrounds = nrounds,
      verbose = verbose
    ), list(...)))
    class(bst) <- c("lgbm", class(bst))
    return(bst)
  }
  stop("Package `lightgbm` not available.")
}

#' @exportS3Method predict lgbm
predict.lgbm <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, data, ...)
}

#' @exportS3Method residuals lgbm
residuals.lgbm <- function(object, response = NULL, data = NULL, ...) {
  preds <- predict(object, data = data, ...)
  .compute_residuals(response, preds)
}
