
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
#' \deqn{E[cox(X, Y | Z)] / \sqrt{E[Var(X | Z)]E[Var(Y | Z)]}},
#' using \code{scale = "cor"}.
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
plm_equiv_test <- function(Y, X, Z, from, to,
                           scale = c("plm", "cov", "cor"), ...) {
  stopifnot(NCOL(X) == 1 && NCOL(Y) == 1)
  stopifnot(from < to)
  scale <- match.arg(scale)

  ### First run standard two-sided test
  tst <- gcm(Y, X, Z, ...)
  n <- NROW(tst$rY)
  ts <- mean(tst$rY * tst$rX)
  vv <- stats::var(tst$rY * tst$rX)
  scale_factor <- switch(
    scale,
    "plm" = c(stats::var(tst$rX)),
    "cov" = 1,
    "cor" = stats::sd(tst$rY) * c(stats::sd(tst$rX))
    )
  t0 <- (from + to) / 2 * scale_factor
  tt <- c("X-squared" = n * (ts - t0)^2 / vv)

  ### Setup
  equiv_margin <- c(from, to) * scale_factor
  width <- diff(equiv_margin)

  ### Compute non-centrality parameter and p-value
  ncp <- (n / vv) * width^2 / 4
  pval <- stats::pchisq(q = tt, df = tst$parameter, ncp = ncp)

  ### Output statistic (chi^2), critical value, p-value
  hyp <- switch(
    scale,
    "plm" = "partial linear effect", # "E[cov(Y, X | Z)] / E[Var(X | Z)]",
    "cov" = "conditional covariance", # "E[cov(Y, X | Z)]",
    "cor" = "conditional correlation" # "E[cov(Y, X | Z)] / sqrt(E[Var(X | Z)] E[Var(Y | Z)])",
  )
  mar <- paste0("in [", from, ", ", to, "]")
  names(mar) <- hyp
  structure(list(
    statistic = tt, p.value = pval, parameter = c(tst$parameter, "ncp" = ncp),
    hypothesis = mar, null.value = mar, alternative = "",
    method = paste0("Partial linear effect equivalence test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    rY = tst$rY, rX = tst$rX), class = c("gcm", "htest"))
}
