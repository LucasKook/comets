#' Equivalence test for the parameter in a partially linear model
#'
#' @details
#' The partially linear model postulates \deqn{Y = X \theta + g(Z) + \epsilon,}
#' and the target of inference is theta. The target is closely related to
#' the conditional covariance between Y and X given Z:
#' \deqn{\theta = E[cov(X, Y | Z)] / E[Var(X | Z)].} The equivalence test (based
#' on the GCM test) tests \eqn{H_0: \theta \not\in [{\tt from}, {\tt to}]} versus
#' \eqn{H_1: \theta \in [{\tt from}, {\tt to}]}. Y, X (and theta) can only be
#' one-dimensional. There are no restrictions on Z. The equivalence test can
#' also be performed on the conditional covariance scale directly (using
#' \code{scale = "cov"}) or on the conditional correlation scale:
#' \deqn{E[cov(X, Y | Z)] / \sqrt{E[Var(X | Z)]E[Var(Y | Z)]}},
#' using \code{scale = "cor"}.
#' The \code{"plm"} and \code{"cor"} tests are experimental and issue a warning.
#' For these scales, the standard error is computed from
#' the influence function of the ratio, accounting for estimation of its
#' denominator. The calculation assumes that the conditional mean regressions
#' and their estimation errors satisfy the conditions needed for asymptotic
#' inference; the score is orthogonal at the true conditional means, but this
#' alone does not guarantee valid inference if one regression is misspecified.
#'
#' @inheritParams gcm
#'
#' @param from Lower bound of the equivalence margin
#' @param to Upper bound of the equivalence margin
#' @param scale Scale on which to specify the equivalence margin. Default
#'     \code{"plm"} corresponds to the partially linear model parameter
#'     described in the details. \code{"cov"} corresponds to the conditional
#'     covariance and \code{"cor"} to conditional correlation which lies in
#'     \eqn{[-1, 1]}.
#' @param ... Further arguments passed to \code{\link{gcm}}
#'
#' @return Object of class '\code{gcm}' and '\code{htest}', including the
#'     estimated parameter in \code{estimate} and its influence-function
#'     standard error in \code{std.error}.
#' @export
#'
#' @examples
#' n <- 150
#' X <- rnorm(n)
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- X^2 + Z[, 2] + rnorm(n)
#' plm_equiv_test(Y, X, Z, from = -1, to = 1)
plm_equiv_test <- function(Y, X, Z, from, to,
                           scale = c("plm", "cov", "cor"), ...) {
  stopifnot(NCOL(X) == 1 && NCOL(Y) == 1)
  stopifnot(length(from) == 1L, length(to) == 1L, is.finite(from), is.finite(to), from < to)
  scale <- match.arg(scale)
  if (scale != "cov") {
    warning("Equivalence testing on the '", scale, "' scale is experimental.", call. = FALSE)
  }

  ### First run standard two-sided test
  tst <- gcm(Y, X, Z, ...)
  rX <- c(tst$rX)
  rY <- c(tst$rY)
  if (length(rX) != length(rY) || any(!is.finite(rX)) || any(!is.finite(rY))) {
    stop("The fitted regressions must produce finite, one-dimensional residuals.", call. = FALSE)
  }
  n <- length(rX)
  xy <- mean(rX * rY)
  if (scale == "cov") {
    estimate <- xy
    influence <- rX * rY - xy
  } else {
    xx <- mean(rX^2)
    yy <- mean(rY^2)
    if (!is.finite(xx) || xx <= 0 || (scale == "cor" && (!is.finite(yy) || yy <= 0))) {
      stop("Residual second moments must be positive and finite.", call. = FALSE)
    }
    if (scale == "plm") {
      estimate <- xy / xx
      influence <- rX * (rY - estimate * rX) / xx
    } else {
      estimate <- xy / sqrt(xx * yy)
      influence <- rX * rY / sqrt(xx * yy) -
        estimate / 2 * (rX^2 / xx + rY^2 / yy)
    }
  }
  vv <- stats::var(influence)
  if (!is.finite(vv) || vv <= 0) stop("Estimated influence-function variance must be positive and finite.", call. = FALSE)
  se <- sqrt(vv / n)
  tt <- c("X-squared" = ((estimate - (from + to) / 2) / se)^2)
  ncp <- ((to - from) / (2 * se))^2
  pval <- unname(stats::pchisq(q = tt, df = 1, ncp = ncp))

  ### Output statistic (chi^2), critical value, p-value
  hyp <- switch(scale,
    "plm" = "partial linear effect", # "E[cov(Y, X | Z)] / E[Var(X | Z)]",
    "cov" = "conditional covariance", # "E[cov(Y, X | Z)]",
    "cor" = "conditional correlation" # "E[cov(Y, X | Z)] / sqrt(E[Var(X | Z)] E[Var(Y | Z)])",
  )
  mar <- paste0("in [", from, ", ", to, "]")
  names(mar) <- hyp
  structure(list(
    statistic = tt, p.value = pval, parameter = c("df" = 1, "ncp" = ncp),
    estimate = stats::setNames(estimate, hyp), std.error = se,
    hypothesis = mar, null.value = mar, alternative = "",
    method = paste0("Partial linear effect equivalence test"),
    data.name = paste0(deparse(match.call(), width.cutoff = 80), collapse = "\n"),
    rY = tst$rY, rX = tst$rX
  ), class = c("gcm", "htest"))
}
