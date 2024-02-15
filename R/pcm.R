#' Projected covariance measure test for conditional mean independence
#'
#' @details
#' The projected covariance measure test tests whether the conditional
#' mean of Y given X and Z depends on X.
#'
#' @references
#' Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2022). The
#' Projected Covariance Measure for assumption-lean variable significance
#' testing. arXiv preprint arXiv:2211.02039. \doi{10.48550/arXiv.2211.02039}
#'
#' @inheritParams gcm
#' @param rep Number of repetitions with which to repeat the PCM test
#' @param est_vhat Logical; whether to estimate the variance functional
#' @param reg Character; regression method that can be one of
#'     \code{"pcm_ranger"} or \code{"pcm_lasso"}.
#' @param mtry Argument passed to \code{ranger}
#' @param ghat_args Arguments passed to \code{reg}
#' @param ... Additional arguments passed to \code{ranger}
#'
#' @importFrom ranger ranger
#'
#' @returns Object of class '\code{pcm}' and '\code{htest}' with the following
#' components:
#' \itemize{
#' \item{\code{statistic}} {The value of the test statistic.}
#' \item{\code{p.value}} {The p-value for the \code{hypothesis}}
#' \item{\code{parameter}} {In case X is multidimensional, this is the degrees of
#'     freedom used for the chi-squared test.}
#' \item{\code{hypothesis}} {Null hypothesis of conditional mean independence.}
#' \item{\code{null.value}} {Null hypothesis of conditional mean independence.}
#' \item{\code{method}} {The string \code{"Projected covariance measure test"}.}
#' \item{\code{data.name}} {A character string giving the name(s) of the data.}
#' \item{\code{check.data}} {A \code{data.frame} containing the residuals for plotting.}
#' }
#'
#' @export
#'
#' @examples
#' n <- 150
#' X <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- rnorm(n) # X[, 2] + Z[, 2] + rnorm(1e3)
#' (pcm1 <- pcm(Y, X, Z))
#' plot(pcm1)
#'
pcm <- function(Y, X, Z, rep = 1, est_vhat = TRUE,
                reg = c("pcm_ranger", "pcm_lasso"),
                ghat_args = NULL, mtry = identity,
                ...) {
  reg <- match.arg(reg)
  if (rep != 1) {
    pcms <- lapply(seq_len(rep), \(iter) {
      pcm(Y = Y, X = X, Z = Z, rep = 1, est_vhat = est_vhat, reg = reg,
          ghat_args = ghat_args, mtry = mtry, ... = ...)
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

  mhats <- predict(mhat, data = Zte)
  dcheck <- data.frame(
    id = setdiff(seq_len(NROW(Y)), idx),
    resid_mhat = (Yte - mhats),
    resid_mhatfhat = (fhats - mhatfhat$predictions)
  )

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
  test <- x$check.data
  if (requireNamespace("ggplot2") && requireNamespace("tidyr") && requireNamespace("ggpubr")) {
    mpl <- \(xx, yy, pdat, ...) {
      ggplot2::ggplot(pdat, ggplot2::aes(y = .data[[yy]], x = .data[[xx]])) +
        ggplot2::geom_point(alpha = 0.3) +
        ggplot2::geom_smooth(se = FALSE, method = "lm") +
        ggplot2::theme_bw()
    }
    p2 <- mpl("resid_mhat", "resid_mhatfhat", test) +
      ggplot2::labs(x = "Residuals f(X, Z) | Z", y = "Residuals Y | Z")
    print(p2)
  }
  return(invisible(p2))
}
