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
#' @param Y Vector or matrix of response values.
#' @param X Matrix or data.frame of covariates.
#' @param Z Matrix or data.frame of covariates.
#' @param alternative A character string specifying the alternative hypothesis,
#'     must be one of \code{"two.sided"} (default), \code{"greater"} or
#'     \code{"less"}. Only applies if \code{type = "quadratic"} and \code{Y} and
#'     \code{X} are one-dimensional.
#' @param reg_YonZ Character string or function specifying the regression for
#'     Y on Z. See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_XonZ Character string or function specifying the regression for
#'     X on Z. See \code{?\link[comets]{regressions}} for more detail.
#' @param args_YonZ A list of named arguments passed to \code{reg_YonZ}.
#' @param args_XonZ A list of named arguments passed to \code{reg_XonZ}.
#' @param type Type of test statistic, either \code{"quadratic"} (default) or
#'     \code{"max"}. If \code{"max"} is specified, the p-value is computed
#'     based on a bootstrap approximation of the null distribution with
#'     \code{B} samples.
#' @param B Number of bootstrap samples. Only applies if \code{type = "max"} is
#'     used.
#' @param coin Logical; whether or not to use the \code{coin} package for
#'     computing the test statistic and p-value. The \code{coin} package
#'     computes variances with n - 1 degrees of freedom.
#'     The default is \code{TRUE}.
#' @param cointrol List; further arguments passed to
#'     \code{\link[coin]{independence_test}}.
#' @param multivariate Character; specifying which regression can handle
#'     multivariate outcomes (\code{"none"}, \code{"YonZ"}, \code{"XonZ"}, or
#'     \code{"both"}). If \code{"none"}, then the regression is run using each
#'     column in Y (or X) as the response.
#' @param ... Additional arguments passed to \code{reg_YonZ}.
#' @param return_fitted_models Logical; whether to return the fitted regressions
#'     (default is \code{FALSE}).
#'
#' @returns Object of class '\code{gcm}' and '\code{htest}' with the following
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
#' @export
#'
#' @examples
#' n <- 1e2
#' X <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(X) <- c("X1", "X2")
#' Z <- matrix(rnorm(2 * n), ncol = 2)
#' colnames(Z) <- c("Z1", "Z2")
#' Y <- X[, 2]^2 + Z[, 2] + rnorm(n)
#' (gcm1 <- gcm(Y, X, Z))
#'
gcm <- function(
    Y, X, Z, alternative = c("two.sided", "less", "greater"),
    reg_YonZ = "rf", reg_XonZ = "rf", args_YonZ = NULL,
    args_XonZ = NULL, type = c("quadratic", "max", "scalar"), B = 499L,
    coin = TRUE, cointrol = list(distribution = "asymptotic"),
    return_fitted_models = FALSE, multivariate = c("none", "YonZ", "XonZ", "both"), ...) {
  Y <- .check_data(Y, "Y")
  X <- .check_data(X, "X")
  Z <- .check_data(Z, "Z")
  alternative <- match.arg(alternative)
  multivariate <- match.arg(multivariate)
  mvYonZ <- multivariate %in% c("YonZ", "both")
  mvXonZ <- multivariate %in% c("XonZ", "both")
  type <- match.arg(type)
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

  if (coin | NCOL(rY) > 1) {
    tst <- do.call("independence_test", c(list(
      rY ~ rX,
      alternative = alternative, teststat = type
    ), cointrol))
    df <- NCOL(rY) * NCOL(rX)
    stat <- coin::statistic(tst)
    pval <- coin::pvalue(tst)
  } else {
    tst <- .gcm(c(rY), rX, alternative = alternative, type = type, B = B)
    df <- tst$df
    stat <- tst$stat
    pval <- tst$pval
  }

  tname <- "X-squared"
  par <- c("df" = df)
  if (type == "max") {
    tname <- "maxT"
    par <- NULL
  } else if (type == "scalar") {
    tname <- "Z"
    par <- NULL
  }
  names(stat) <- tname

  models <- if (return_fitted_models) {
    list(reg_YonZ = mY, reg_XonZ = mX)
  } else {
    NULL
  }

  structure(list(
    statistic = stat, p.value = pval, parameter = par,
    hypothesis = c("E[cov(Y, X | Z)]" = "0"),
    null.value = c("E[cov(Y, X | Z)]" = "0"), alternative = alternative,
    method = paste0("Generalized covariance measure test"),
    data.name = deparse(match.call(), width.cutoff = 80),
    rY = rY, rX = rX, models = models
  ), class = c("gcm", "htest"))
}

# Helpers -----------------------------------------------------------------

.multi_regression <- function(Y, X, reg, args, rfm, mv) {
  if (mv) {
    m <- do.call(reg, c(list(y = Y, x = X), args))
    r <- stats::residuals(m, response = Y, data = X)
  } else {
    res <- apply(Y, 2, \(tY) {
      m <- do.call(reg, c(list(y = tY, x = X), args))
      r <- stats::residuals(m, response = tY, data = X)
      if (rfm) list(m = m, r = r) else r
    }, simplify = FALSE)
    if (rfm) {
      r <- do.call("cbind", lapply(res, \(x) x[["r"]]))
      m <- lapply(res, \(x) x[["m"]])
    } else {
      r <- do.call("cbind", res)
      m <- NULL
    }
  }
  return(list(models = m, residuals = r))
}

.check_data <- function(x, mode = c("Y", "X", "Z"), test = "gcm") {
  mode <- match.arg(mode)
  if (mode == "Y") {
    N <- NROW(x)
    if (!is.matrix(x) & !is.data.frame(x)) {
      ret <- c(x)
      if (is.factor(ret) && length(levels(ret)) > 2 && test == "pcm") {
        stop("Only binary factors are allowed for Y.")
      }
      return(ret)
    } else {
      if (NCOL(x) > 1 && test == "pcm") {
        stop("Please provide Y as a vector.")
      }
      .check_data(x, mode = "X")
    }
  }
  if ("tibble" %in% class(x)) {
    x <- as.data.frame(x)
  }
  if (NCOL(x) == 1) {
    x <- as.matrix(x)
  }
  if (is.null(colnames(x))) {
    colnames(x) <- paste0(mode, seq_len(NCOL(x)))
  }
  x
}

#' @importFrom stats pnorm
.gcm <- function(
    rY, rX, alternative = "two.sided", type = "quadratic", B = 499L) {
  dY <- NCOL(rY)
  dX <- NCOL(rX)
  nn <- NROW(rY)
  RR <- rY * rX
  if (dY > 1 || dX > 1) {
    if (type == "quadratic") {
      sigma <- crossprod(RR) / nn - tcrossprod(colMeans(RR))
      eig <- eigen(sigma)
      if (min(eig$values) < .Machine$double.eps) {
        warning("`vcov` of test statistic is not invertible")
      }
      siginvhalf <- eig$vectors %*% diag(eig$values^(-1 / 2)) %*%
        t(eig$vectors)
      tstat <- siginvhalf %*% colSums(RR) / sqrt(nn)
      stat <- sum(tstat^2)
      pval <- stats::pchisq(stat, df = dX, lower.tail = FALSE)
    } else {
      tRR <- t(RR)
      mRR <- rowMeans(tRR)
      tRR <- tRR / sqrt((rowMeans(tRR^2) - mRR^2))
      stat <- max(abs(mRR)) * sqrt(nn)
      sim <- apply(abs(tRR %*% matrix(
        stats::rnorm(nn * B), nn, B
      )), 2, max) / sqrt(nn)
      pval <- (sum(sim >= stat) + 1) / (B + 1)
    }
  } else {
    R.sq <- RR^2
    meanR <- mean(RR)
    stat <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
    pval <- switch(alternative,
      "two.sided" = 2 * stats::pnorm(-abs(stat)),
      "greater" = 1 - stats::pnorm(stat),
      "less" = stats::pnorm(stat)
    )
    stat <- stat^2
  }
  list("stat" = stat, "pval" = pval, df = dX)
}

.rm_int <- function(x) {
  if (all(x[, 1] == 1)) {
    return(x[, -1L, drop = FALSE])
  }
  x
}

#' @importFrom stats terms
.get_terms <- function(formula) {
  if (is.null(formula)) {
    return(NULL)
  }
  atms <- stats::terms(formula)
  tms <- attr(atms, "term.labels")
  resp <- all.vars(formula)[1]
  ridx <- grep("|", tms, fixed = TRUE)
  tms[ridx] <- paste0("(", tms[ridx], ")")
  ie <- grep(":", tms, value = TRUE)
  me <- grep(":", tms, value = TRUE, invert = TRUE)
  list(
    all = tms, me = me, ie = ie, response = resp, terms = atms,
    fml = formula
  )
}

# Diagnostics -------------------------------------------------------------

#' Plotting methods for COMETs
#'
#' @rdname plot.comet
#'
#' @param x Object of class '\code{gcm}', '\code{pcm}', or '\code{wgcm}'.
#' @param plot Logical; whether to print the plot (default: \code{TRUE}).
#' @param ... Currently ignored.
#'
#' @exportS3Method plot gcm
plot.gcm <- function(x, plot = TRUE, ...) {
  if (requireNamespace("ggplot2") & requireNamespace("tidyr") & requireNamespace("dplyr")) {
    .data <- NULL
    pd <- tidyr::pivot_longer(data.frame(rY = unname(x$rY), rX = unname(x$rX)),
      dplyr::starts_with("rX"),
      names_to = "nX",
      values_to = "rX"
    )
    if (NCOL(x$rY > 1)) {
      pd <- tidyr::pivot_longer(pd, dplyr::starts_with("rY"),
        names_to = "nY", values_to = "rY"
      )
    } else {
      pd$nY <- "rY.1"
    }
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(
      x = .data[["rX"]], y = .data[["rY"]],
      color = interaction(.data[["nY"]], .data[["nX"]]),
      linetype = .data[["nY"]]
    )) +
      ggplot2::geom_point(alpha = 0.3, show.legend = FALSE) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Residuals X | Z", y = "Residuals Y | Z")
    if (plot) print(p1)
    return(invisible(p1))
  }
  stop("Package `ggplot2` not available.")
}

.mm <- function(preds, data) {
  .rm_int(stats::model.matrix(stats::reformulate(preds), data = data))
}

#' GCM test with pre-computed residuals
#'
#' @param rY Vector or matrix of response values.
#' @param rX Matrix or data.frame of covariates.
#' @param type Type of test statistic, either \code{"quadratic"} (default) or
#'     \code{"max"}. If \code{"max"} is specified, the p-value is computed
#'     based on a bootstrap approximation of the null distribution with
#'     \code{B} samples.
#' @param alternative A character string specifying the alternative hypothesis,
#'     must be one of \code{"two.sided"} (default), \code{"greater"} or
#'     \code{"less"}. Only applies if \code{type = "quadratic"} and \code{Y} and
#'     \code{X} are one-dimensional.
#' @param ... Further arguments passed to \code{\link[coin]{independence_test}()}.
#'
#' @returns Object of class '\code{gcm}' and '\code{htest}' with the following
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
#' @export
rgcm <- function(
    rY, rX, alternative = "two.sided",
    type = c("quadratic", "max", "scalar"), ...) {
  type <- match.arg(type)
  tst <- coin::independence_test(rY ~ rX,
    alternative = alternative, teststat = type, ...
  )
  df <- NULL
  switch(type,
    "scalar" = {
      stat <- c("Z" = coin::statistic(tst))
    },
    "quadratic" = {
      stat <- c("X-squared" = coin::statistic(tst))
      df <- tst@statistic@df
    },
    "max" = {
      stat <- c("maxT" = coin::statistic(tst))
    }
  )
  structure(list(
    statistic = stat, p.value = coin::pvalue(tst), parameter = df,
    hypothesis = c("E[cov(Y, X | Z)]" = "0"),
    null.value = c("E[cov(Y, X | Z)]" = "0"), alternative = alternative,
    method = paste0("Generalized covariance measure test"),
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    rY = rY, rX = rX
  ), class = c("gcm", "htest"))
}
