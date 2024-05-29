
#' Equivalence test for the parameter in a partially linear model
#'
#' @details
#' The partially linear model postulates \deqn{Y = X \theta + g(Z) + \epsilon,}
#' and the target of inference is theta. The target is closely related to
#' the conditional covariance between Y and X given Z:
#' \deqn{\theta = E[cov(X, Y | Z)] / E[Var(X | Z)].} The equivalence test (based
#' on the GCM test) tests \eqn{H_0: \theta \not\in [{\tt from}, {\tt to}]} versus
#' \eqn{H_1: \theta \in [{\tt from}, {\tt to}]}. Y, X (and theta) can only be
#' one-dimensional. There are no restrictions on Z.
#'
#' @inheritParams gcm
#'
#' @param from Lower bound of the equivalence margin
#' @param to Upper bound of the equivalence margin
#' @param ... Further arguments passed to \code{\link{gcm}}
#'
#' @return Object of class '\code{gcm}' and '\code{htest}'
#' @export
#'
#' @examples
#' n <- 150
#' X <- rnorm(n)
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- X^2 + Z[, 2] + rnorm(n)
#' plm_equiv_test(Y, X, Z, from = -1, to = 1)
plm_equiv_test <- function(Y, X, Z, from, to, ...) {
  stopifnot(NCOL(X) == 1 && NCOL(Y) == 1)
  stopifnot(from < to)

  ### First run standard two-sided test
  tst <- gcm(Y, X, Z, ...)
  n <- NROW(tst$rY)
  ts <- mean(tst$rY * tst$rX)
  vv <- stats::var(tst$rY * tst$rX)
  t0 <- (from + to) / 2 * stats::var(tst$rX)
  tt <- c("X-squared" = n * (ts - t0)^2 / vv)

  ### Setup
  equiv_margin <- c(from, to) * stats::var(c(tst$rX))
  width <- diff(equiv_margin)

  ### Function to compute the critical value for the equivalence test
  cval <- \(s, level, w = width) {
    stats::qchisq(p = level, df = 1, ncp = (s * w / 2)^2)
  }

  ### Compute p-value via root-finding
  pval <- stats::uniroot(\(alp) cval(sqrt(n) / sqrt(vv), level = alp) - tt,
                         interval = c(0, 1))$root

  ### Output statistic (chi^2), critical value, p-value
  structure(list(
    statistic = tt, p.value = pval, parameter = tst$parameter,
    hypothesis = c("E[cov(Y, X | Z)] / E[Var(X | Z)]" = paste0(
      "in [", from, ", ", to, "]")),
    null.value = c("E[cov(Y, X | Z)] / E[Var(X | Z)]" = paste0(
      "in [", from, ", ", to, "]")),
    alternative = "",
    method = paste0("Partially linear effect equivalence test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    rY = tst$rY, rX = tst$rX), class = c("gcm", "htest"))
}
