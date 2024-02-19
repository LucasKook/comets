#' Covariance measure tests with formula interface
#'
#' @details
#' Formula-based interface for the generalised and projected covariance measure
#' tests.
#'
#' @param formula Formula of the form \code{Y ~ X | Z} for testing Y independent
#'     of X given Z.
#' @param data Data.frame containing the variables in \code{formula}.
#' @param test Character string; \code{"gcm"} or \code{"pcm"}.
#' @param ... Additional arguments passed to \code{test}.
#'
#' @return Object of class \code{"gcm"} or \code{"pcm"} and \code{"htest"}.
#'     See \code{\link{gcm}} and \code{\link{pcm}} for details.
#' @export
#'
#'
#' @examples
#' tn <- 3e2
#' df <- data.frame(y = rnorm(tn), x1 = rnorm(tn), x2 = rnorm(tn), z = rnorm(tn))
#' comet(y ~ x1 + x2 | z, data = df, test = "gcm")
comet <- function(formula, data, test = c("gcm", "pcm"), ...) {
  fm <- Formula::as.Formula(formula)
  Y <- stats::model.response(model.frame(fm, data))
  X <- .rm_int(stats::model.matrix(fm, data, rhs = 1))
  Z <- .rm_int(stats::model.matrix(fm, data, rhs = 2))
  test <- match.fun(match.arg(test))
  test(Y = Y, X = X, Z = Z, ...)
}
