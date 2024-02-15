#' Generalised covariance measure test using random forests
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
#'     \code{"less"}
#' @param ... Additional arguments passed to ranger
#'
#' @returns Object of class '\code{gcm}' and '\code{htest}' with the following
#' components:
#' \itemize{
#' \item{\code{statistic}} {The value of the test statistic.}
#' \item{\code{p.value}} {The p-value for the \code{hypothesis}}
#' \item{\code{parameter}} {In case X is multidimensional, this is the degrees of
#'     freedom used for the chi-squared test.}
#' \item{\code{hypothesis}} {String specifying the null hypothesis .}
#' \item{\code{null.value}} {String specifying the null hypothesis.}
#' \item{\code{method}} {The string \code{"Generalised covariance measure test"}.}
#' \item{\code{data.name}} {A character string giving the name(s) of the data.}
#' \item{\code{rY}} {Residuals for the Y on Z regression.}
#' \item{\code{rX}} {Residuals for the X on Z regression.}
#' }
#' @export
#'
#' @examples
#' X <- matrix(rnorm(2e3), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2e3), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- rnorm(1e3) # X[, 2] + Z[, 2] + rnorm(1e3)
#' (gcm1 <- gcm(Y, X, Z))
#'
gcm <- function(Y, X, Z, alternative = c("two.sided", "less", "greater"), ...) {
  alternative <- match.arg(alternative)
  args <- if (length(list(...)) > 0) list(...) else NULL
  YZ <- do.call("pcm_ranger", c(list(y = Y, x = Z), args))
  XZ <- apply(as.data.frame(X), 2, \(tX) {
    do.call("pcm_ranger", c(list(y = tX, x = Z), args))
  })
  rY <- Y - predict(YZ, data = Z)
  preds <- lapply(XZ, predict.pcm_ranger, data = Z)
  rX <- X - do.call("cbind", preds)
  stat <- .gcm(rY, rX)
  pval <- .compute_normal_pval(stat, alternative)

  if (!is.null(df <- attr(pval, "df"))) {
    tname <- "X-squared"
    par <- c("df" = df)
  } else {
    tname <- "Z"
    par <- NULL
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

.gcm <- function (r, e) {
  dR <- NCOL(r)
  dE <- NCOL(e)
  nn <- NROW(r)
  if (dR > 1 || dE > 1) {
    R_mat <- matrix(r, nrow = nn, ncol = dE) * e
    sigma <- crossprod(R_mat)/nn - tcrossprod(colMeans(R_mat))
    eig <- eigen(sigma)
    if (min(eig$values) < .Machine$double.eps)
      warning("`vcov` of test statistic is not invertible")
    siginvhalf <- eig$vectors %*% diag(eig$values^(-1/2)) %*%
      t(eig$vectors)
    tstat <- siginvhalf %*% colSums(R_mat)/sqrt(nn)
    stat <- structure(sum(tstat^2), df = dR * dE)
  }
  else {
    R <- r * e
    R.sq <- R^2
    meanR <- mean(R)
    stat <- sqrt(nn) * meanR/sqrt(mean(R.sq) - meanR^2)
  }
  stat
}

#' @importFrom stats pnorm
.compute_normal_pval <- function(stat, alternative) {
  if (!is.null(df <- attr(stat, "df")))
    return(stats::pchisq(stat, df = df, lower.tail = FALSE))
  switch(
    alternative,
    "two.sided" = 2 * stats::pnorm(-abs(stat)),
    "greater" = stats::pnorm(-abs(stat)),
    "less" = stats::pnorm(abs(stat))
  )
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