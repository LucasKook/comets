#' Weighted Generalised covariance measure test
#'
#' @details
#' The weighted generalised covariance measure test tests whether a weighted
#' version of the conditional covariance of Y and X given Z is zero.
#'
#' @references
#' Scheidegger, C., Hörrmann, J., & Bühlmann, P. (2022). The weighted
#' generalised covariance measure. Journal of Machine Learning Research,
#' 23(273), 1-68.
#'
#' @inherit gcm
#'
#' @param Y Vector of response values. Can be supplied as a numeric vector or
#'     a single column matrix.
#' @param reg_wfun Character string or function specifying the regression for
#'     estimating the weighting function.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param args_wfun Additional arguments passed to \code{reg_XonZ}.
#' @param frac Relative size of train split.
#' @param ... Additional arguments passed to \code{reg_YonZ}.
#'
#' @returns Object of class '\code{wgcm}' and '\code{htest}' with the following
#' components:
#' \item{\code{statistic}}{The value of the test statistic.}
#' \item{\code{p.value}}{The p-value for the \code{hypothesis}}
#' \item{\code{parameter}}{In case X is multidimensional, this is the degrees of
#'     freedom used for the chi-squared test.}
#' \item{\code{hypothesis}}{String specifying the null hypothesis .}
#' \item{\code{null.value}}{String specifying the null hypothesis.}
#' \item{\code{method}}{The string \code{"Generalised covariance measure test"}.}
#' \item{\code{data.name}}{A character string giving the name(s) of the data.}
#' \item{\code{rY}}{Residuals for the Y on Z regression.}
#' \item{\code{rX}}{Weighted residuals for the X on Z regression.}
#'
#' @export
#'
#' @examples
#' n <- 150
#' X <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- X[, 2]^2 + Z[, 2] + rnorm(n)
#' (wgcm1 <- wgcm(Y, X, Z))
#'
wgcm <- function(Y, X, Z, reg_YonZ = "rf", reg_XonZ = "rf", reg_wfun = "rf",
                 args_XonZ = NULL, args_wfun = NULL, frac = 0.5,
                 B = 499L, coin = FALSE, cointrol = NULL, ...) {
  Y <- .check_data(Y, "Y")
  X <- .check_data(X, "X")
  Z <- .check_data(Z, "Z")
  alternative <- "greater"
  args <- if (length(list(...)) > 0) list(...) else NULL

  ### Sample splitting
  dsp <- .split_sample(Y, X, Z, frac = frac)
  Ytr <- dsp$Ytr
  Xtr <- dsp$Xtr
  Ztr <- dsp$Ztr
  Yte <- dsp$Yte
  Xte <- dsp$Xte
  Zte <- dsp$Zte

  ### Estimate weight function and compute weights
  wYZ <- do.call(reg_YonZ, c(list(y = Ytr, x = Ztr), args))
  wrY <- stats::residuals(wYZ, response = Ytr, data = Ztr)
  wrX <- apply(as.data.frame(Xtr), 2, \(tX) {
    mX <- do.call(reg_XonZ, c(list(y = tX, x = Ztr), args_XonZ))
    stats::residuals(mX, response = tX, data = Ztr)
  })
  RP <- c(wrY) * wrX
  W <- apply(as.data.frame(RP), 2, \(tRP) {
    mZ <- do.call(reg_wfun, c(list(y = tRP, x = Ztr), args_wfun))
    sign(predict(mZ, data = Zte))
  })

  ### GCM on test data with weights
  YZ <- do.call(reg_YonZ, c(list(y = Yte, x = Zte), args))
  rY <- stats::residuals(YZ, response = Yte, data = Zte)
  rX <- apply(as.data.frame(Xte), 2, \(tX) {
    mX <- do.call(reg_XonZ, c(list(y = tX, x = Zte), args_XonZ))
    stats::residuals(mX, response = tX, data = Zte)
  })
  if (coin) {
    tst <- do.call("independence_test", c(list(
      rY ~ rX, alternative = "greater", teststat = "max",
      distribution = coin::approximate(B))))
    df <- NCOL(rY) * NCOL(rX)
    stat <- coin::statistic(tst)
    pval <- coin::pvalue(tst)
  } else {
    tst <- .gcm(rY, as.matrix(rX * W), alternative = "greater",
                type = "max", B = B)
    df <- tst$df
    stat <- tst$stat
    pval <- tst$pval
  }

  tname <- ifelse(df == 1, "Z", "|Z|")
  par <- NULL
  names(stat) <- tname

  structure(list(
    statistic = stat, p.value = pval, parameter = par,
    hypothesis = c("E[w(Z) cov(Y, X | Z)]" = "0"),
    null.value = c("E[w(Z) cov(Y, X | Z)]" = "0"), alternative = alternative,
    method = paste0("Weighted generalized covariance measure test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    rY = rY, rX = rX * W), class = c("wgcm", "htest"))

}

# Vis ---------------------------------------------------------------------

#' @exportS3Method plot wgcm
plot.wgcm <- function(x, ...) {
  .data <- NULL
  pd <- tidyr::pivot_longer(data.frame(rY = x$rY, rX = unname(x$rX)),
                            dplyr::starts_with("rX"))
  if (requireNamespace("ggplot2")) {
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data[["value"]] ,
                                           y = .data[["rY"]],
                                           color = .data[["name"]])) +
      ggplot2::geom_point(alpha = 0.3, show.legend = FALSE) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Weighted residuals X | Z", y = "Residuals Y | Z")
    print(p1)
  }
  return(invisible(p1))
}
