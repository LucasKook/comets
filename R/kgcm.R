#' Kernel generalised covariance measure test
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
#'     multivariate outcomes (\code{"none"}, or \code{"XonZ"}).
#'     If \code{"none"}, then the regression is run using each column in X as
#'     the response.
#' @param return_fitted_models Logical; whether to return the fitted regressions
#'     (default is \code{FALSE}).
#' @param bandwidth Numeric; value of the bandwidth for the Gaussian kernel.
#'     Defaults to \code{NULL}, corresponding to the median heuristic.
#'     For constant Z the default produces a constant kernel; if the median
#'     distance is zero for nonconstant Z, provide a positive bandwidth.
#' @param return_kernel Logical; retain the Gram matrix in the result for
#'     plotting (default \code{TRUE}). Setting it to \code{FALSE} reduces the
#'     size of the returned object, but not peak memory during computation.
#' @param batch_size Maximum number of bootstrap draws multiplied by the Gram
#'     matrix at once. Larger batches may be faster but require more memory.
#' @param ... Additional arguments passed to \code{reg_YonZ}.
#'
#' @returns Object of class '\code{kgcm}' and '\code{htest}' with the following
#' components:
#' \item{\code{statistic}}{The value of the test statistic.}
#' \item{\code{p.value}}{The p-value for the \code{hypothesis}}
#' \item{\code{hypothesis}}{String specifying the null hypothesis.}
#' \item{\code{null.value}}{String specifying the null hypothesis.}
#' \item{\code{method}}{The string \code{"Kernel generalised covariance measure test"}.}
#' \item{\code{data.name}}{A character string giving the name(s) of the data.}
#' \item{\code{rY}}{Residuals for the Y on Z regression.}
#' \item{\code{rX}}{Residuals for the X on Z regression.}
#' \item{\code{K}}{Gram matrix of the Gaussian kernel evaluated at Z.}
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
  multivariate = c("none", "XonZ"),
  bandwidth = NULL, return_kernel = TRUE, batch_size = 32L, ...
) {
  Y <- .check_data(Y, "Y", "pcm")
  X <- .check_data(X, "X")
  Z <- .check_data(Z, "Z")
  .validate_observations(Y = Y, X = X, Z = Z)
  .validate_bootstrap(B)
  if (length(batch_size) != 1L || !is.numeric(batch_size) || !is.finite(batch_size) || batch_size < 1 || batch_size != floor(batch_size)) {
    stop("`batch_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(return_kernel) || length(return_kernel) != 1L || is.na(return_kernel)) {
    stop("`return_kernel` must be TRUE or FALSE.", call. = FALSE)
  }
  multivariate <- match.arg(multivariate)
  mvXonZ <- multivariate %in% c("XonZ", "both")
  args <- if (length(list(...)) > 0) list(...) else NULL
  args <- c(args_YonZ, args)
  mY <- do.call(reg_YonZ, c(list(y = Y, x = Z), args))
  rY <- stats::residuals(mY, response = Y, data = Z)
  XZ <- .multi_regression(X, Z, reg_XonZ, args_XonZ, return_fitted_models, mvXonZ)
  rX <- XZ[["residuals"]]
  mX <- XZ[["models"]]

  Zm <- as.matrix(Z)
  if (is.null(bandwidth)) {
    h <- medianheur(Zm, NROW(Z), NCOL(Z))
    if (h == 0 && all(sweep(Zm, 2L, Zm[1L, ]) == 0)) h <- 1
  } else if (is.numeric(bandwidth) && length(bandwidth) == 1L) {
    h <- bandwidth
  } else {
    stop("Supply a valid bandwidth (`NULL` for the median heuristic, or a numeric value).")
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("Bandwidth must be positive and finite; supply `bandwidth` if the median distance is zero.", call. = FALSE)
  }
  K <- gaussian_gram(Zm, h, NROW(Z), NCOL(Z))
  RP <- c(rY) * rX
  obs <- 1 / NROW(Z) * sum(RP * (K %*% RP))
  smpl <- numeric(B)
  for (start in seq.int(1L, B, by = batch_size)) {
    count <- min(batch_size, B - start + 1L)
    multipliers <- matrix(sample(c(-1, 1), NROW(Z) * count, TRUE), NROW(Z), count)
    values <- do.call(cbind, lapply(seq_len(count), \(i) RP * multipliers[, i]))
    quadratic <- colSums(values * (K %*% values)) / NROW(Z)
    smpl[seq.int(start, length.out = count)] <- colSums(matrix(quadratic, nrow = NCOL(RP)))
  }
  pval <- (1 + sum(smpl >= obs)) / (1 + B)

  models <- if (return_fitted_models) {
    list(reg_YonZ = mY, reg_XonZ = mX)
  } else {
    NULL
  }

  structure(list(
    statistic = c("Psi" = obs), p.value = pval, parameter = NULL,
    hypothesis = stats::setNames("0", .hypothesis(rY, "w(Z) ")),
    null.value = stats::setNames("0", .hypothesis(rY, "w(Z) ")), alternative = "two.sided",
    method = paste0("Kernel generalised covariance measure test"),
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    rY = rY, rX = rX, K = if (return_kernel) K else NULL, models = models
  ), class = c("kgcm", "htest"))
}

# Vis ---------------------------------------------------------------------

#' @rdname plot.comet
#' @exportS3Method plot kgcm
plot.kgcm <- function(x, plot = TRUE, ...) {
  .data <- NULL
  if (is.null(x$K)) stop("The Gram matrix was not retained; set `return_kernel = TRUE` when fitting.", call. = FALSE)
  if (requireNamespace("ggplot2") & requireNamespace("tidyr") & requireNamespace("dplyr")) {
    pd <- tidyr::pivot_longer(
      data.frame(rY = x$rY, rX = unname(x$K %*% x$rX)),
      dplyr::starts_with("rX")
    )
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
    return(invisible(p1))
  }
  stop("Package `ggplot2` not available.")
}
