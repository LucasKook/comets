#' Implemented regression methods
#' @rdname regressions
#' @param y Vector (or matrix) of response values.
#' @param x Design matrix of predictors.
#' @param ... Additional arguments passed to the underlying regression method.
#'    In case of \code{"rf"}, \code{"tuned_rf"}, \code{"survforest"} and
#'    \code{"qrf"}, this is \code{\link[ranger]{ranger}}. In case of
#'    \code{"lasso"} and \code{"ridge"}, this is \code{\link[glmnet]{glmnet}}.
#'    In case of \code{"cox"}, this is \code{\link[survival]{coxph}}. In case
#'    of \code{"xgb"} and \code{"tuned_xgb"} these are the parameters passed
#'    to \code{\link[xgboost]{xgb.train}}. In case of \code{"lgbm"}, this is
#'    \code{lgb.train} from the \code{lightgbm} package.
#' @details
#' The implemented choices are \code{"rf"} for random forests as implemented in
#' ranger, \code{"lasso"} for cross-validated Lasso regression (using the
#' penalty specified in \code{s}), \code{"ridge"}
#' for cross-validated ridge regression (using the penalty specified in \code{s}),
#' \code{"cox"} for the Cox proportional
#' hazards model as implemented in survival, \code{"qrf"} or \code{"survforest"}
#' for quantile and survival random forests, respectively. The option
#' \code{"postlasso"} option refers to a cross-validated LASSO (using the
#' penalty specified in \code{s}) and subsequent OLS regression. The \code{"lrm"}
#' option implements a standard linear regression model. The \code{"xgb"} and
#' \code{"tuned_xgb"} options require the \code{xgboost} package and the
#' \code{"lgbm"} option requires the \code{lightgbm} package.
#'
#' The \code{"tuned_rf"} regression method tunes the \code{mtry} and
#' \code{max.depth} parameters in \code{\link[ranger]{ranger}} out-of-bag.
#' The \code{"tuned_xgb"} regression method uses k-fold cross-validation to
#' tune the \code{nrounds}, \code{eta} and \code{max_depth} parameters in
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
#' For ordered factors with more than two levels, \code{"rf"} and
#' \code{"tuned_rf"} estimate the conditional distribution of the response
#' with a probability forest. The residuals are then computed from the
#' estimated probabilities \eqn{p_k(Z) = P(Y = k \mid Z)} and cumulative
#' probabilities \eqn{P_k(Z) = P(Y \le k \mid Z)}, \eqn{k = 1, \dots, K},
#' according to \code{ordinal_residuals}:
#' \code{"probability"} (default) gives the \eqn{K - 1}-dimensional residuals
#' \eqn{1(Y = k) - p_k(Z)}, \eqn{k = 2, \dots, K}, which have power against
#' general alternatives; \code{"mean"} gives the one-dimensional residual
#' \eqn{Y - \sum_k k p_k(Z)} using the integer scores \eqn{1, \dots, K},
#' which has power against location shifts; \code{"score"} gives the
#' one-dimensional score residual of a proportional odds model, i.e., the
#' score with respect to \eqn{\gamma} in
#' \eqn{\mathrm{logit} P(Y \le k \mid Z) = \mathrm{logit} P_k(Z) - \gamma}
#' at \eqn{\gamma = 0}, which has the closed form
#' \eqn{P_{Y-1}(Z) + P_Y(Z) - 1 = P(Y^\prime < Y \mid Z) -
#' P(Y^\prime > Y \mid Z)} and has power against proportional odds
#' alternatives.
#'
#' @param ordinal_residuals Type of residual computed for ordered factor
#'     responses with more than two levels, one of \code{"probability"}
#'     (default), \code{"mean"} or \code{"score"}. See Details. Ignored for
#'     other response types.
#'
rf <- function(y, x, ordinal_residuals = c("probability", "mean", "score"),
               ...) {
  ordinal_residuals <- match.arg(ordinal_residuals)
  args <- list(...)
  y <- .drop_unused_levels(y)
  if (length(unique(y)) == 2) {
    y <- factor(y)
  }
  if (is.factor(y)) {
    args$probability <- TRUE
  }
  rf <- do.call("ranger", c(list(y = y, x = x), args))
  rf$ordinal_residuals <- ordinal_residuals
  rf$response_levels <- levels(y)
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
  response <- .align_levels(response, object$response_levels)
  type <- object$ordinal_residuals
  if (is.ordered(response) && nlevels(response) > 2 && !is.null(type) &&
    type != "probability") {
    class(object) <- class(object)[-1]
    probs <- predict(object, data = data)$predictions
    return(.ordinal_residuals(response, probs, type))
  }
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

### Internal residuals function for ordered factors, `probs` is the n x K
### matrix of estimated conditional probabilities P(Y = k | Z)
.ordinal_residuals <- function(y, probs, type = c("mean", "score")) {
  type <- match.arg(type)
  K <- nlevels(y)
  if (NCOL(probs) != K) {
    stop("Number of predicted probabilities does not match number of levels.")
  }
  k <- as.integer(y)
  if (type == "mean") {
    return(k - c(probs %*% seq_len(K)))
  }
  ### Cumulative probabilities P_0 = 0, P_1, ..., P_K = 1
  cprobs <- cbind(0, probs %*% upper.tri(diag(K), diag = TRUE))
  idx <- seq_along(k)
  cprobs[cbind(idx, k + 1)] + cprobs[cbind(idx, k)] - 1
}

### Drop unused levels of a factor response with a warning
.drop_unused_levels <- function(y) {
  if (!is.factor(y)) {
    return(y)
  }
  unused <- setdiff(levels(y), levels(droplevels(y)))
  if (length(unused) > 0) {
    warning(
      "Dropping unused level(s) of the response: ",
      paste(unused, collapse = ", "), call. = FALSE
    )
    y <- droplevels(y)
  }
  y
}

### Restrict levels of a factor response to those seen during training
.align_levels <- function(y, lvls) {
  if (!is.factor(y) || is.null(lvls) || identical(levels(y), lvls)) {
    return(y)
  }
  if (!all(y %in% lvls)) {
    stop(
      "Response contains level(s) not seen during training: ",
      paste(setdiff(unique(as.character(y)), lvls), collapse = ", ")
    )
  }
  factor(y, levels = lvls, ordered = is.ordered(y))
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
  pred <- \(i) sum(K[i, ] * as.numeric(object$response <= object$response[i]))
  ipreds <- sapply(seq_along(object$response), pred)
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
  s = "lambda.min", ...
) {
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
                     ordinal_residuals = c("probability", "mean", "score"),
                     ...) {
  ordinal_residuals <- match.arg(ordinal_residuals)
  args <- list(...)
  y <- .drop_unused_levels(y)
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
  rf$ordinal_residuals <- ordinal_residuals
  rf$response_levels <- levels(y)
  class(rf) <- c("rf", class(rf))
  rf
}

# Boosting ---------------------------------------------------------------

#' @rdname regressions
#' @param nrounds See \code{\link[xgboost]{xgb.train}}.
#' @param verbose See \code{\link[xgboost]{xgb.train}}.
xgb <- function(y, x, nrounds = 2L, verbose = 0L, ...) {
  if (requireNamespace("xgboost")) {
    bst <- xgboost::xgb.train(
      params = list(...), data = xgboost::xgb.DMatrix(x, label = y),
      nrounds = nrounds, verbose = verbose
    )
    class(bst) <- c("xgb", class(bst))
    return(bst)
  }
  stop("Package `xgboost` not available.")
}

#' @rdname regressions
#' @param etas Values for \code{eta} to cross-validate. See
#'     \code{\link[xgboost]{xgb.train}}.
#' @param max_depths Values for \code{max_depth} to cross-validate. See
#'     \code{\link[xgboost]{xgb.train}}.
#' @param metrics See \code{\link[xgboost]{xgb.cv}}.
#' @param nfold Number of folds for \code{nfold}-cross validation.
#' @param folds Specify folds for cross validation.
tuned_xgb <- function(y, x, nfold, folds, etas = c(0.1, 0.5, 1),
                      max_depths = 1:5, nrounds = c(2, 10, 50),
                      verbose = 0, metrics = list("rmse"), ...) {
  if (requireNamespace("xgboost")) {
    args <- list(
      data = xgboost::xgb.DMatrix(x, label = y), verbose = verbose,
      metrics = metrics
    )
    if (missing(folds) && missing(nfold)) {
      nfold <- 5
    }
    if (!missing(nfold)) {
      args$nfold <- nfold
    }
    if (!missing(folds)) {
      args$folds <- folds
    }
    cvres <- lapply(etas, \(teta) {
      lapply(max_depths, \(tmd) {
        lapply(nrounds, \(tnr) {
          args$nrounds <- tnr
          args$params <- c(list(eta = teta, max_depth = tmd), list(...))
          cv <- do.call(xgboost::xgb.cv, args)
          ### Take min over boosting iterations
          err <- min(cv$evaluation_log[[paste0("test_", metrics[[1]], "_mean")]])
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
    bst <- do.call(lightgbm::lgb.train, c(list(
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
