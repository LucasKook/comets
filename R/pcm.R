#' Conditional mean independence test
#'
#' @param Y Numeric; response
#' @param X Numeric; covariates
#' @param Z Numceric; covariates
#' @param rep Number of repetitions
#' @param est_vhat Estimate variance functional
#' @param ... Additional arguments passed to \code{reg}
#' @param reg Character; regression method
#' @param mtry Argument passed to \code{ranger}
#' @param ghat_args Arguments passed to reg
#' @param do.check Save check data
#'
#' @importFrom ranger ranger
#'
#' @return Object of class '\code{htest}'
#' @export
#'
#' @examples
#' X <- matrix(rnorm(2e3), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2e3), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- rnorm(1e3) # X[, 2] + Z[, 2] + rnorm(1e3)
#' (pcm1 <- pcm(Y, X, Z, ghat_args = list(mtry = NULL, max.depth = NULL),
#'     est_vhat = TRUE, do.check = TRUE))
#' plot(pcm1)
#'
pcm <- function(Y, X, Z, rep = 1, est_vhat = TRUE,
                reg = c("pcm_ranger", "pcm_lasso"),
                ghat_args = NULL, mtry = identity,
                do.check = FALSE, ...) {
  reg <- match.arg(reg)
  if (rep != 1) {
    pcms <- lapply(seq_len(rep), \(iter) {
      pcm(Y = Y, X = X, Z = Z, rep = 1, est_vhat = est_vhat, reg = reg,
          ghat_args = ghat_args, mtry = mtry, do.check = FALSE, ... = ...)
    })
    stat <- mean(unlist(lapply(pcms, \(tst) tst$statistic)))
    pval <- pnorm(stat, lower.tail = FALSE)
    return(structure(list(
      statistic = c("Z" = stat), p.value = pval,
      hypothesis = c("E[Y | X, Z]" = "E[Y | Z]"),
      null.value = c("E[Y | X, Z]" = "E[Y | Z]"), alternative = "two.sided",
      method = paste0("Projected covariance measure test (K = ", rep, " repetitions)"),
      all_tests = pcms, data.name = deparse(match.call(), width.cutoff = 80)),
      class = c("pcm", "htest")))
  }
  ### Sample splitting
  idx <- sample.int(NROW(Y), ceiling(NROW(Y) / 2))
  ### Split 1
  Ytr <- Y[idx]
  Xtr <- data.frame(X)[idx, , drop = FALSE]
  Ztr <- data.frame(Z)[idx, , drop = FALSE]
  ### Split 2
  Yte <- Y[-idx]
  Xte <- data.frame(X)[-idx, , drop = FALSE]
  Zte <- data.frame(Z)[-idx, , drop = FALSE]

  ### Obtain hat{h}
  ghat <- do.call(reg, c(list(y = Ytr, x = cbind(Xtr, Ztr)), ghat_args))
  mtilde <- ranger(x = Ztr, y = pghat <- predict(ghat, data = cbind(Xtr, Ztr)),
                   mtry = mtry, ...)
  htilde <- \(X, Z) {
    predict(ghat, data = cbind(X, Z)) -
      predict(mtilde, data = Z)$predictions
  }
  rho <- mean((Ytr - mtilde$predictions) * predict(ghat, data = cbind(Xtr, Ztr)))
  hhat <- \(X, Z) sign(rho) * htilde(X, Z)

  ### Obtain hat{v}
  if (est_vhat) {
    vtilde <- ranger(x = cbind(Xtr, Ztr), y = (sqr <- (Ytr - predict(
      ghat, data = cbind(Xtr, Ztr)))^2), mtry = mtry, ...)
    a <- function(c) mean(sqr / (pmax(vtilde$predictions, 0) + c))
    chat <- if (a(0) < 1) 0 else stats::uniroot(\(c) a(c) - 1, c(0, 10), extendInt = "yes")$root
    vhat <- \(X, Z) pmax(predict(vtilde, data = cbind(X, Z))$predictions, 0) + chat
  }
  else
    vhat <- \(X, Z) 1

  ### Obtain residuals for test
  fhat <- \(X, Z) hhat(X, Z) / vhat(X, Z)
  mhatfhat <- ranger(x = Zte, y = (fhats <- fhat(Xte, Zte)), mtry = mtry, ...)
  mhat <- do.call(reg, c(list(y = Yte, x = Zte), ghat_args))

  ### Test
  L <- (Yte - predict(mhat, data = Zte)) * (fhats - mhatfhat$predictions)
  stat <- sqrt(length(idx)) * mean(L) / sqrt(mean(L^2) - mean(L)^2)
  if (is.nan(stat)) stat <- -Inf
  pval <- pnorm(stat, lower.tail = FALSE)

  dcheck <- if (do.check) {
    pmtilde <- predict(mtilde, data = Ztr)$predictions
    # pvtilde <- predict(vtilde, data = cbind(Xtr, Ztr))$predictions
    pghat_test <- predict(ghat, data = cbind(Xte, Zte))
    mhats <- predict(mhat, data = Zte)
    list(
      train = data.frame(
        id = idx,
        fitted_ghat = pghat,
        resid_ghat = Ytr - pghat,
        fitted_mtilde = pmtilde,
        resid_mtilde = pghat - pmtilde,
        # fitted_vtilde = pvtilde,
        # resid_vtilde = sqr - pvtilde,
        hhat = hhat(Xtr, Ztr),
        vhat = vhat(Xtr, Ztr),
        MSE_YZ = (Ytr - predict(mhat, data = Ztr))^2,
        MSE_YXZ = (Ytr - pghat)^2
      ),
      test = data.frame(
        id = setdiff(seq_len(NROW(Y)), idx),
        fitted_ghat = pghat_test,
        resid_ghat = Yte - pghat_test,
        hhat = hhat(Xte, Zte),
        vhat = vhat(Xte, Zte),
        fhat = fhats,
        resid_mhat = (Yte - mhats),
        mhat = mhats,
        resid_mhatfhat = (fhats - mhatfhat$predictions),
        mhatfhat = mhats,
        MSE_YXZ = (Yte - pghat_test)^2,
        MSE_YZ = (Yte - mhats)^2
      )
    )
  } else NULL

  structure(list(
    statistic = c("Z" = stat), p.value = pval,
    hypothesis = c("E[Y | X, Z]" = "E[Y | Z]"),
    null.value = c("E[Y | X, Z]" = "E[Y | Z]"), alternative = "two.sided",
    method = paste0("Projected covariance measure test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    check.data = dcheck), class = c("pcm", "htest"))
}

# Regressions -------------------------------------------------------------

pcm_ranger <- function(y, x, ...) {
  args <- list(...)
  if (length(unique(y)) == 2) {
    y <- as.factor(y)
    args$probability <- TRUE
  }
  rf <- do.call("ranger", c(list(y = y, x = x), args))
  class(rf) <- c("pcm_ranger", class(rf))
  rf
}

predict.pcm_ranger <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  preds <- predict(object, data = data)$predictions
  if (object$treetype == "Probability estimation")
    preds <- preds[, 2]
  preds
}

#' @importFrom glmnet cv.glmnet
pcm_lasso <- function(y, x, ...) {
  obj <- cv.glmnet(y = y, x = as.matrix(x), ...)
  class(obj) <- c("pcm_lasso", class(obj))
  obj
}

predict.pcm_lasso <- function(object, data = NULL, ...) {
  class(object) <- class(object)[-1]
  predict(object, newx = as.matrix(data), s = object$lambda.min)[, 1]
}

# Diagnostics -------------------------------------------------------------

#' @exportS3Method plot pcm
plot.pcm <- function(x, ...) {
  .data <- NULL
  train <- x$check.data$train
  test <- x$check.data$test
  mse_train <- data.frame(MSE_YZ = train$MSE_YZ, MSE_YXZ = train$MSE_YXZ, set = "train")
  mse_test <- data.frame(MSE_YZ = test$MSE_YZ, MSE_YXZ = test$MSE_YXZ, set = "test")
  mses <- tidyr::pivot_longer(dplyr::bind_rows(mse_train, mse_test),
                              dplyr::starts_with("MSE"))
  if (is.null(train))
    return("Nothing to plot. Consider running `pcm` with `do.check = TRUE`")
  if (requireNamespace("ggplot2") && requireNamespace("tidyr") && requireNamespace("ggpubr")) {
    mpl <- \(xx, yy, pdat, ...) {
      ggplot2::ggplot(pdat, ggplot2::aes(y = .data[[yy]], x = .data[[xx]])) +
        ggplot2::geom_point(alpha = 0.3) +
        ggplot2::geom_smooth(se = FALSE, method = "lm") +
        ggplot2::theme_bw()
    }
    # ptrain <- lapply(c("ghat", "mtilde", "vtilde"), \(tx) {
    #   mpl(paste0("fitted_", tx), paste0("resid_", tx), train, "direction")
    # })
    # ptest <- apply(data.frame(y = c("resid_mhatfhat", "resid_mhat"),
    #                           x = c("mhatfhat", "resid_mhatfhat")), 1,
    #                \(x) mpl(x[2], x[1], test), simplify = FALSE)
    # p1 <- mpl("mhatfhat", "resid_mhatfhat", test) +
    #   ggplot2::labs(x = "f(X,Z)", y = "Residuals f(X,Z) | Z")
    p2 <- mpl("resid_mhat", "resid_mhatfhat", test) +
      ggplot2::labs(x = "Residuals f(X, Z) | Z", y = "Residuals Y | Z")
    p3 <- ggplot2::ggplot(mses, ggplot2::aes(y = .data[["value"]],
                                             x = .data[["name"]],
                                             color = .data[["set"]])) +
      ggplot2::geom_violin(width = 0.5, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::geom_boxplot(width = 0.3, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::geom_jitter(alpha = 0.1, position = ggplot2::position_dodge(width = 0.7)) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_log10() +
      ggplot2::labs(y = "MSE contributions", x = "Regression")
    print(p2)
  }
  return(invisible(p2))
}
