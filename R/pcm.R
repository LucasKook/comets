#' Projected covariance measure test for conditional mean independence
#'
#' @details
#' The projected covariance measure test tests whether the conditional
#' mean of Y given X and Z is independent of X.
#'
#' @references
#' Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2024). The
#' Projected Covariance Measure for assumption-lean variable significance
#' testing. The Annals of Statistics, 52(6), 2851-2878.
#' \doi{10.1214/24-AOS2447}
#'
#' @inheritParams gcm
#' @param Y Vector of response values. Can be supplied as a numeric vector or
#'     a single column matrix.
#' @param rep Number of repetitions with which to repeat the PCM test
#' @param est_vhat Logical; whether to estimate the variance functional
#' @param reg_YonXZ Character string or function specifying the regression
#'     for Y on X and Z, default is \code{"rf"} for random forest.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_YonZ Character string or function specifying the regression
#'     for Y on Z, default is \code{"rf"} for random forest.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_YhatonZ Character string or function specifying the regression
#'     for the predicted values of \code{reg_YonXZ} on Z, default is \code{"rf"}
#'     for random forest.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_VonXZ Character string or function specifying the regression
#'     for estimating the conditional variance of Y given X and Z, default
#'     is \code{"rf"} for random forest.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param reg_RonZ Character string or function specifying the regression
#'     for the estimated transformation of Y, X, and Z on Z, default is
#'     \code{"rf"} for random forest.
#'     See \code{?\link[comets]{regressions}} for more detail.
#' @param args_YonXZ A list of named arguments passed to \code{reg_YonXZ}.
#' @param args_YonZ A list of named arguments passed to \code{reg_YonZ}.
#' @param args_YhatonZ A list of named arguments passed to \code{reg_YhatonZ}.
#' @param args_VonXZ A list of named arguments passed to \code{reg_VonXZ}.
#' @param args_RonZ A list of named arguments passed to \code{reg_RonZ}.
#' @param frac Relative size of train split.
#' @param indices A numeric vector of indices specifying the observations used
#'     for estimating the estimating the direction (the other observations will
#'     be used for computing the final test statistic). Default is \code{NULL}
#'     and the indices will be generated randomly using \code{frac}.
#'     When using \code{rep} larger than 1, a list (of length \code{rep}) of
#'     indices can be supplied.
#' @param return_fitted_models Logical; whether to return the fitted regressions
#'     (default is \code{FALSE}).
#' @param ... Additional arguments currently ignored.
#'
#' @importFrom ranger ranger
#' @importFrom coin independence_test
#'
#' @returns Object of class '\code{pcm}' and '\code{htest}' with the following
#' components:
#' \item{\code{statistic}}{The value of the test statistic.}
#' \item{\code{p.value}}{The p-value for the \code{hypothesis}}
#' \item{\code{hypothesis}}{Null hypothesis of conditional mean independence.}
#' \item{\code{null.value}}{Null hypothesis of conditional mean independence.}
#' \item{\code{method}}{The string \code{"Projected covariance measure test"}.}
#' \item{\code{data.name}}{A character string giving the name(s) of the data.}
#' \item{\code{check.data}}{A \code{data.frame} containing the residuals for plotting.}
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
#' (pcm1 <- pcm(Y, X, Z))
#'
pcm <- function(
  Y, X, Z, rep = 1, est_vhat = TRUE, reg_YonXZ = "rf",
  reg_YonZ = "rf", reg_YhatonZ = "rf", reg_VonXZ = "rf",
  reg_RonZ = "rf", args_YonXZ = NULL, args_YonZ = NULL,
  args_YhatonZ = NULL, args_VonXZ = NULL, args_RonZ = NULL,
  frac = 0.5, indices = NULL, coin = FALSE, cointrol = NULL,
  return_fitted_models = FALSE,
  ...
) {
  call <- match.call()
  ### Data checks
  Y <- .check_data(Y, "Y", "pcm")
  X <- .check_data(X, "X", "pcm")
  Z <- .check_data(Z, "Z", "pcm")
  n <- .validate_observations(Y = Y, X = X, Z = Z)
  if (length(rep) != 1L || !is.numeric(rep) || !is.finite(rep) || rep < 1 || rep != floor(rep)) {
    stop("`rep` must be a positive integer.", call. = FALSE)
  }
  if (rep > 1L && !is.null(indices) && (!is.list(indices) || length(indices) != rep)) {
    stop("Please supply a list of indices of length `rep`.", call. = FALSE)
  }
  if (rep == 1L) .validate_split(n, frac, indices)
  if (rep > 1L) {
    for (i in seq_len(rep)) .validate_split(n, frac, if (is.null(indices)) NULL else indices[[i]])
  }
  ### Multiple splits
  if (rep != 1) {
    if (!is.null(indices) && length(indices) != rep) {
      stop("Please supply a list of indices of length `rep`.")
    }
    pcms <- lapply(seq_len(rep), \(iter) {
      pcm(
        Y = Y, X = X, Z = Z, rep = 1, est_vhat = est_vhat,
        reg_YonXZ = reg_YonXZ, reg_YonZ = reg_YonZ,
        reg_YhatonZ = reg_YhatonZ, reg_VonXZ = reg_VonXZ,
        reg_RonZ = reg_RonZ, args_YonXZ = args_YonXZ,
        args_YonZ = args_YonZ, args_YhatonZ = args_YhatonZ,
        args_VonXZ = args_VonXZ, args_RonZ = args_RonZ,
        frac = frac, indices = indices[[iter]], coin = coin,
        cointrol = cointrol,
        return_fitted_models = return_fitted_models, ...
      )
    })
    stat <- mean(unlist(lapply(pcms, \(tst) tst$statistic)))
    pval <- pnorm(stat, lower.tail = FALSE)
    dcheck <- do.call("rbind", lapply(1:rep, \(x) {
      data.frame(iter = x, pcms[[x]]$check.data)
    }))
  } else {
    ### Sample splitting
    dsp <- .split_sample(Y, X, Z, frac = frac, indices = indices)
    Ytr <- dsp$Ytr
    Xtr <- dsp$Xtr
    Ztr <- dsp$Ztr
    Yte <- dsp$Yte
    Xte <- dsp$Xte
    Zte <- dsp$Zte
    idx <- dsp$idx

    ### Obtain hat{h}
    ghat <- do.call(reg_YonXZ, c(list(y = Ytr, x = cbind(Xtr, Ztr)), args_YonXZ))
    pghat <- predict(ghat, data = cbind(Xtr, Ztr))
    mtilde <- do.call(reg_YhatonZ, c(list(x = Ztr, y = pghat), args_YhatonZ))
    htilde <- \(X, Z) predict(ghat, data = cbind(X, Z)) - predict(mtilde, data = Z)
    rho <- mean(stats::residuals(mtilde, response = Ytr, data = Ztr) * pghat)
    hhat <- \(X, Z) sign(rho) * htilde(X, Z)

    ### Obtain hat{v}
    if (est_vhat) {
      sqr <- stats::residuals(ghat, response = Ytr, data = cbind(Xtr, Ztr))^2
      vtilde <- do.call(reg_VonXZ, c(list(x = cbind(Xtr, Ztr), y = sqr), args_VonXZ))
      vtrain <- pmax(predict(vtilde, data = cbind(Xtr, Ztr)), 0)
      if (any(!is.finite(vtrain)) || any(!is.finite(sqr))) {
        stop("Variance regression produced non-finite values.", call. = FALSE)
      }
      a <- function(c) {
        denom <- vtrain + c
        ratio <- sqr / denom
        ratio[denom == 0 & sqr == 0] <- 0
        mean(ratio)
      }
      chat <- if (a(0) < 1) 0 else stats::uniroot(\(c) a(c) - 1, c(0, 10), extendInt = "yes")$root
      vhat <- function(X, Z) {
        v <- pmax(predict(vtilde, data = cbind(X, Z)), 0) + chat
        if (any(!is.finite(v) | v <= 0)) stop("Estimated test variance must be positive and finite.", call. = FALSE)
        v
      }
    } else {
      vtilde <- NULL
      vhat <- \(X, Z) 1
    }

    ### Obtain residuals for test
    fhat <- \(X, Z) hhat(X, Z) / vhat(X, Z)
    fhats <- fhat(Xte, Zte)
    mhatfhat <- do.call(reg_RonZ, c(list(x = Zte, y = fhats), args_RonZ))
    mhat <- do.call(reg_YonZ, c(list(y = Yte, x = Zte), args_YonZ))

    ### Test
    rY <- stats::residuals(mhat, response = Yte, data = Zte)
    rT <- fhats - predict(mhatfhat, data = Zte)

    if (coin) {
      tst <- do.call("independence_test", c(list(
        rY ~ rT,
        alternative = "greater", teststat = "scalar"
      ), cointrol))
      stat <- coin::statistic(tst)
      pval <- coin::pvalue(tst)
    } else {
      L <- rY * rT
      if (any(!is.finite(L))) stop("PCM residual products must be finite.", call. = FALSE)
      variance <- mean((L - mean(L))^2)
      if (!is.finite(variance) || variance <= 0) stop("PCM residual products must have positive finite variance.", call. = FALSE)
      stat <- sqrt(NROW(Yte)) * mean(L) / sqrt(variance)
      pval <- pnorm(stat, lower.tail = FALSE)
    }

    dcheck <- data.frame(
      id = setdiff(seq_len(NROW(Y)), idx),
      rY = rY, rT = rT, iter = 1
    )
  }

  models <- if (return_fitted_models && rep == 1) {
    list(
      reg_YonXZ = ghat, reg_YhatonZ = mtilde, reg_VonXZ = vtilde,
      reg_RonZ = mhatfhat, reg_YonZ = mhat
    )
  } else if (return_fitted_models && rep > 1) {
    lapply(pcms, \(x) x[["models"]])
  } else {
    NULL
  }

  structure(list(
    statistic = c("Z" = stat), p.value = pval,
    hypothesis = c("E[Y | X, Z]" = "E[Y | Z]"),
    null.value = c("E[Y | X, Z]" = "E[Y | Z]"), alternative = "two.sided",
    method = paste0("Projected covariance measure test"),
    data.name = paste0(deparse(match.call()), collapse = "\n"),
    check.data = dcheck, rep = rep, models = models
  ), class = c("pcm", "htest"))
}

# Helpers -----------------------------------------------------------------

.validate_split <- function(n, frac, indices = NULL) {
  if (is.null(indices)) {
    if (length(frac) != 1L || !is.numeric(frac) || !is.finite(frac) || frac <= 0 || frac >= 1) {
      stop("`frac` must be strictly between 0 and 1.", call. = FALSE)
    }
    size <- ceiling(frac * n)
  } else {
    if (!is.numeric(indices) || anyNA(indices) || any(!is.finite(indices)) ||
      any(indices != floor(indices)) || any(indices < 1 | indices > n) || anyDuplicated(indices)) {
      stop("`indices` must contain unique, valid integer row indices.", call. = FALSE)
    }
    size <- length(indices)
  }
  if (size < 2L || n - size < 2L) stop("Each split must contain at least two observations.", call. = FALSE)
  invisible(NULL)
}

.split_sample <- function(Y, X, Z, frac = 0.5, indices = NULL) {
  idx <- indices
  if (is.null(idx)) {
    idx <- sample.int(NROW(Y), ceiling(frac * NROW(Y)))
  }
  ### Split 1
  Ytr <- Y[idx]
  Xtr <- as.matrix(data.frame(X)[idx, , drop = FALSE])
  Ztr <- as.matrix(data.frame(Z)[idx, , drop = FALSE])
  ### Split 2
  Yte <- Y[-idx]
  Xte <- as.matrix(data.frame(X)[-idx, , drop = FALSE])
  Zte <- as.matrix(data.frame(Z)[-idx, , drop = FALSE])
  list(Ytr = Ytr, Xtr = Xtr, Ztr = Ztr, Yte = Yte, Xte = Xte, Zte = Zte, idx = idx)
}

# Diagnostics -------------------------------------------------------------

#' @rdname plot.comet
#' @exportS3Method plot pcm
plot.pcm <- function(x, plot = TRUE, ...) {
  .data <- NULL
  test <- x$check.data
  if (requireNamespace("ggplot2")) {
    mpl <- \(xx, yy, pdat, ...) {
      ggplot2::ggplot(pdat, ggplot2::aes(
        y = .data[[yy]], x = .data[[xx]],
        color = factor(.data[["iter"]])
      )) +
        ggplot2::geom_point(alpha = 0.3, show.legend = FALSE) +
        ggplot2::geom_smooth(se = FALSE, method = "lm", show.legend = FALSE) +
        ggplot2::theme_bw()
    }
    p2 <- mpl("rT", "rY", test) +
      ggplot2::labs(x = "Residuals f(X, Z) | Z", y = "Residuals Y | Z")
    if (plot) print(p2)
    return(invisible(p2))
  }
  stop("Package `ggplot2` not available.")
}
