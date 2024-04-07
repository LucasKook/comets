#' Generalised covariance measure test
#'
#' @details
#' The generalised covariance measure test tests whether the conditional
#' covariance of Y and X given Z is zero.
#'
#' @references
#' Rajen D. Shah, Jonas Peters "The hardness of conditional independence testing
#' and the generalised covariance measure," The Annals of Statistics, 48(3),
#' 1514-1538. \doi{10.1214/19-aos1857}
#'
#' @param Y Vector of response values. Can be supplied as a numeric vector or
#'     a single column matrix.
#' @param X Matrix or data.frame of covariates.
#' @param Z Matrix or data.frame of covariates.
#' @param alternative A character string specifying the alternative hypothesis,
#'     must be one of \code{"two.sided"} (default), \code{"greater"} or
#'     \code{"less"}. Only applies if \code{type = "quadratic"} and \code{X}
#'     is one-dimensional.
#' @param reg_YonZ Character string or function specifying the regression for
#'     Y on Z.
#' @param reg_XonZ Character string or function specifying the regression for
#'     X on Z.
#' @param args_XonZ Additional arguments passed to \code{reg_XonZ}.
#' @param type Type of test statistic, either \code{"quadratic"} (default) or
#'     \code{"max"}. If \code{"max"} is specified, the p-value is computed
#'     based on a bootstrap approximation of the null distribution with
#'     \code{B} samples.
#' @param B Number of bootstrap samples. Only applies if \code{type = "max"} is
#'     used.
#' @param ... Additional arguments passed to \code{reg_YonZ}.
#'
#' @returns Object of class '\code{gcm}' and '\code{htest}' with the following
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
#' \item{\code{rX}}{Residuals for the X on Z regression.}
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
#' (gcm1 <- gcm(Y, X, Z))
#'
gcm <- function(Y, X, Z, alternative = c("two.sided", "less", "greater"),
                reg_YonZ = "rf", reg_XonZ = "rf", args_XonZ = NULL,
                type = c("quadratic", "max"), B = 499L, ...) {
  Y <- .check_data(Y, "Y")
  X <- .check_data(X, "X")
  Z <- .check_data(Z, "Z")
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  args <- if (length(list(...)) > 0) list(...) else NULL
  YZ <- do.call(reg_YonZ, c(list(y = Y, x = Z), args))
  rY <- stats::residuals(YZ, response = Y, data = Z)
  rX <- apply(as.data.frame(X), 2, \(tX) {
    mX <- do.call(reg_XonZ, c(list(y = tX, x = Z), args_XonZ))
    stats::residuals(mX, response = tX, data = Z)
  })
  tst <- .gcm(rY, rX, alternative = alternative, type = type, B = B)
  df <- tst$df
  stat <- tst$stat
  pval <- tst$pval

  tname <- "Z"
  par <- NULL
  if (type == "quadratic" && tst$df > 1) {
    tname <- "X-squared"
    par <- c("df" = df)
  } else if (type == "max") {
    tname <- "|Z|"
  }
  names(stat) <- tname

  structure(list(
    statistic = stat, p.value = pval, parameter = par,
    hypothesis = c("E[cov(Y, X | Z)]" = "0"),
    null.value = c("E[cov(Y, X | Z)]" = "0"), alternative = alternative,
    method = paste0("Generalized covariance measure test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    rY = rY, rX = rX), class = c("gcm", "htest"))

}

# Helpers -----------------------------------------------------------------

.compute_residuals <- function(y, pred) {
  if (is.factor(y) && length(levels(y)) == 2)
    y <- as.numeric(y) - 1
  y - pred
}

.check_data <- function(x, mode = c("Y", "X", "Z")) {
  mode <- match.arg(mode)
  if (mode == "Y") {
    N <- NROW(x)
    ret <- c(x)
    if (is.factor(ret) && length(levels(ret)) > 2)
      stop("Only binary factors are allowed for Y.")
    if (N != NROW(ret))
      stop("Please provide Y as a vector.")
    return(ret)
  }
  if ("tibble" %in% class(x))
    x <- as.data.frame(x)
  if (NCOL(x) == 1)
    x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste0(mode, seq_len(NCOL(x)))
  x
}

#' @importFrom stats pnorm
.gcm <- function(
    rY, rX, alternative = "two.sided", type = c("quadratic", "max"), B = 499L
) {
  dY <- NCOL(rY)
  dX <- NCOL(rX)
  nn <- NROW(rY)
  RR <- rY * rX
  if (dY > 1 || dX > 1) {
    if (type == "quadratic") {
      sigma <- crossprod(RR)/nn - tcrossprod(colMeans(RR))
      eig <- eigen(sigma)
      if (min(eig$values) < .Machine$double.eps)
        warning("`vcov` of test statistic is not invertible")
      siginvhalf <- eig$vectors %*% diag(eig$values^(-1/2)) %*%
        t(eig$vectors)
      tstat <- siginvhalf %*% colSums(RR) / sqrt(nn)
      stat <- sum(tstat^2)
      pval <- stats::pchisq(stat, df = dX, lower.tail = FALSE)
    } else {
      tRR <- t(RR)
      mRR <- rowMeans(tRR)
      tRR <- tRR/sqrt((rowMeans(tRR^2) - mRR^2))
      stat <- max(abs(mRR)) * sqrt(nn)
      sim <- apply(abs(tRR %*% matrix(
        stats::rnorm(nn * B), nn, B)), 2, max) / sqrt(nn)
      pval <- (sum(sim >= stat) + 1)/(B + 1)
    }
  }
  else {
    R.sq <- RR^2
    meanR <- mean(RR)
    stat <- sqrt(nn) * meanR/sqrt(mean(R.sq) - meanR^2)
    pval <- switch(
      alternative,
      "two.sided" = 2 * stats::pnorm(-abs(stat)),
      "greater" = 1 - stats::pnorm(stat),
      "less" = stats::pnorm(stat)
    )
  }
  list("stat" = stat, "pval" = pval, df = dX)
}

.rm_int <- function(x) {
  if (all(x[, 1] == 1))
    return(x[, -1L, drop = FALSE])
  x
}

#' @importFrom stats terms
.get_terms <- function(formula) {
  if (is.null(formula))
    return(NULL)
  atms <- stats::terms(formula)
  tms <- attr(atms, "term.labels")
  resp <- all.vars(formula)[1]
  ridx <- grep("|", tms, fixed = TRUE)
  tms[ridx] <- paste0("(", tms[ridx], ")")
  ie <- grep(":", tms, value = TRUE)
  me <- grep(":", tms, value = TRUE, invert = TRUE)
  list(all = tms, me = me, ie = ie, response = resp, terms = atms,
       fml = formula)
}

# Ranger ------------------------------------------------------------------

#' @importFrom stats model.response model.frame
.ranger <- function(formula, data, ...) {
  response <- stats::model.response(stats::model.frame(formula, data))
  is_factor <- is.factor(response)
  tms <- .get_terms(formula)
  resp <- if (is_factor)
    .rm_int(stats::model.matrix(~ response, contrasts.arg = list(
      "response" = "contr.treatment")))
  else response
  tmp <- list(data = data, response = resp, is_factor = is_factor)
  if (identical(tms$me, character(0))) {
    if (is_factor)
      return(structure(c(list(mean = base::colMeans(resp)), tmp),
                       class = "ranger"))
    else return(structure(c(list(mean = mean(as.numeric(response))),
                            tmp), class = "ranger"))
  }
  ret <- ranger::ranger(formula, data, probability = is_factor, ...)
  structure(c(ret, tmp), class = "ranger")
}

#' @importFrom stats predict
residuals.ranger <- function(object, newdata = NULL, newy = NULL, ...) {
  if (is.null(newdata))
    newdata <- object$data
  if (!is.null(newy))
    newy <- if (object$is_factor)
      .rm_int(stats::model.matrix(~ newy, contrasts.arg = list(
        "newy" = "contr.treatment")))
  else newy
  if (is.null(newy))
    newy <- object$response
  if (!is.null(object$mean))
    return(newy - object$mean)
  preds <- stats::predict(object, data = newdata)$predictions
  if (object$is_factor)
    preds <- preds[, -1]
  unname(newy - preds)
}

# Diagnostics -------------------------------------------------------------

#' @exportS3Method plot gcm
plot.gcm <- function(x, ...) {
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
      ggplot2::labs(x = "Residuals X | Z", y = "Residuals Y | Z")
    print(p1)
  }
  return(invisible(p1))
}

.mm <- function(preds, data)
  .rm_int(stats::model.matrix(stats::reformulate(preds), data = data))
