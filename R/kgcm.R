#' Kernelized generalised covariance measure test
#'
#' @details
#' The kernelized generalised covariance measure test tests whether the weighted
#' conditional covariance of Y and X given Z is zero.
#'
#' @references
#' Fernández, T., & Rivera, N. (2024). A general framework for the analysis of
#' kernel-based tests. Journal of Machine Learning Research, 25(95), 1-40.
#'
#' @param Y Vector of response values.
#' @param X Matrix or data.frame of covariates.
#' @param Z Matrix or data.frame of covariates.
#' @param reg_YonZ Character string or function specifying the regression for
#'     Y on Z. See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_XonZ Character string or function specifying the regression for
#'     X on Z. See \code{?\link[comets]{regressions}} for more detail.
#' @param args_YonZ A list of named arguments passed to \code{reg_YonZ}.
#' @param args_XonZ A list of named arguments passed to \code{reg_XonZ}.
#' @param B Number of wild bootstrap samples.
#' @param multivariate Character; specifying which regression can handle
#'     multivariate outcomes (\code{"none"}, \code{"YonZ"}, \code{"XonZ"}, or
#'     \code{"both"}). If \code{"none"}, then the regression is run using each
#'     column in Y (or X) as the response.
#' @param return_fitted_models Logical; whether to return the fitted regressions
#'     (default is \code{FALSE}).
#' @param bandwidth Numeric; value of the bandwidth for the Gaussian kernel.
#'     Defaults to \code{NULL}, corresponding to the median heuristic.
#' @param ... Currently ignored
#'
#' @returns Object of class '\code{kgcm}' and '\code{htest}' with the following
#' components:
#' \item{\code{statistic}}{The value of the test statistic.}
#' \item{\code{p.value}}{The p-value for the \code{hypothesis}}
#' \item{\code{parameter}}{In case X is multidimensional, this is the degrees of
#'     freedom used for the chi-squared test.}
#' \item{\code{hypothesis}}{String specifying the null hypothesis.}
#' \item{\code{null.value}}{String specifying the null hypothesis.}
#' \item{\code{method}}{The string \code{"Generalised covariance measure test"}.}
#' \item{\code{data.name}}{A character string giving the name(s) of the data.}
#' \item{\code{rY}}{Residuals for the Y on Z regression.}
#' \item{\code{rX}}{Residuals for the X on Z regression.}
#' \item{\code{models}}{List of fitted regressions if \code{return_fitted_models} is \code{TRUE}.}
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib comets, .registration = TRUE
#' @export
#'
#' @examples
#' n <- 1e2
#' X <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- X[, 2]^2 + Z[, 2] + rnorm(n)
#' (gcm1 <- kgcm(Y, X, Z))
#'
kgcm <- function(
    Y, X, Z, reg_YonZ = "rf", reg_XonZ = "rf", args_YonZ = NULL,
    args_XonZ = NULL, B = 499L, return_fitted_models = FALSE,
    multivariate = c("none", "YonZ", "XonZ", "both"),
    bandwidth = NULL, ...) {
  Y <- .check_data(Y, "Y", "pcm")
  X <- .check_data(X, "X")
  Z <- .check_data(Z, "Z")
  multivariate <- match.arg(multivariate)
  mvYonZ <- multivariate %in% c("YonZ", "both")
  mvXonZ <- multivariate %in% c("XonZ", "both")
  args <- if (length(list(...)) > 0) list(...) else NULL
  args <- c(args_YonZ, args)
  if ("matrix" %in% class(Y)) {
    YZ <- .multi_regression(Y, Z, reg_YonZ, args, return_fitted_models, mvYonZ)
    rY <- YZ[["residuals"]]
    mY <- YZ[["models"]]
  } else {
    mY <- do.call(reg_YonZ, c(list(y = Y, x = Z), args))
    rY <- stats::residuals(mY, response = Y, data = Z)
  }
  XZ <- .multi_regression(X, Z, reg_XonZ, args_XonZ, return_fitted_models, mvXonZ)
  rX <- XZ[["residuals"]]
  mX <- XZ[["models"]]

  if (is.null(bandwidth)) {
    h <- medianheur(Z, NROW(Z), NCOL(Z))
  } else if (is.numeric(bandwidth)) {
    h <- bandwidth
  } else {
    stop("Supply a valid bandwidth (`NULL` for the median heuristic, or a numeric value).")
  }
  K <- gaussian_gram(Z, h, NROW(Z), NCOL(Z))
  RP <- rY * rX
  obs <- 1 / NROW(Z) * c(t(RP) %*% K %*% RP)
  smpl <- sapply(1:B, \(b) {
    RM <- sample(c(-1, 1), NROW(Z), TRUE)
    RPM <- RP * RM
    1 / NROW(Z) * c(t(RPM) %*% K %*% RPM)
  })
  pval <- (1 + sum(smpl > obs)) / (1 + B)

  models <- if (return_fitted_models) {
    list(reg_YonZ = mY, reg_XonZ = mX)
  } else {
    NULL
  }

  structure(list(
    statistic = c("Psi" = obs), p.value = pval, parameter = NULL,
    hypothesis = c("E[w(Z) cov(Y, X | Z)]" = "0"),
    null.value = c("E[w(Z) cov(Y, X | Z)]" = "0"), alternative = "two.sided",
    method = paste0("Kernelized generalized covariance measure test"),
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    rY = rY, rX = rX, K = K, models = models
  ), class = c("kgcm", "htest"))
}

# Vis ---------------------------------------------------------------------

#' @rdname plot.comet
#' @exportS3Method plot kgcm
plot.kgcm <- function(x, plot = TRUE, ...) {
  .data <- NULL
  pd <- tidyr::pivot_longer(
    data.frame(rY = x$rY, rX = unname(x$K %*% x$rX)),
    dplyr::starts_with("rX")
  )
  if (requireNamespace("ggplot2")) {
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(
      x = .data[["value"]],
      y = .data[["rY"]],
      color = .data[["name"]]
    )) +
      ggplot2::geom_point(alpha = 0.3, show.legend = FALSE) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Weighted residuals X | Z", y = "Residuals Y | Z")
    if (plot) print(p1)
  }
  return(invisible(p1))
}
