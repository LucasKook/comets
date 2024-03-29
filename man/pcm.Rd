% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcm.R
\name{pcm}
\alias{pcm}
\title{Projected covariance measure test for conditional mean independence}
\usage{
pcm(
  Y,
  X,
  Z,
  rep = 1,
  est_vhat = TRUE,
  reg_YonXZ = "rf",
  reg_YonZ = "rf",
  args_YonXZ = NULL,
  args_YonZ = NULL,
  mtry = identity,
  ...
)
}
\arguments{
\item{Y}{Vector of response values. Can be supplied as a numeric vector or
a single column matrix.}

\item{X}{Matrix or data.frame of covariates.}

\item{Z}{Matrix or data.frame of covariates.}

\item{rep}{Number of repetitions with which to repeat the PCM test}

\item{est_vhat}{Logical; whether to estimate the variance functional}

\item{reg_YonXZ}{Character string or function specifying the regression
for Y on X and Z, default is \code{"rf"} for random forest.}

\item{reg_YonZ}{Character string or function specifying the regression
for Y on Z, default is \code{"rf"} for random forest.}

\item{args_YonXZ}{Arguments passed to \code{reg}}

\item{args_YonZ}{Arguments passed to \code{reg}}

\item{mtry}{Argument passed to \code{ranger}}

\item{...}{Additional arguments passed to \code{ranger}}
}
\value{
Object of class '\code{pcm}' and '\code{htest}' with the following
components:
\item{\code{statistic}}{The value of the test statistic.}
\item{\code{p.value}}{The p-value for the \code{hypothesis}}
\item{\code{parameter}}{In case X is multidimensional, this is the degrees of
    freedom used for the chi-squared test.}
\item{\code{hypothesis}}{Null hypothesis of conditional mean independence.}
\item{\code{null.value}}{Null hypothesis of conditional mean independence.}
\item{\code{method}}{The string \code{"Projected covariance measure test"}.}
\item{\code{data.name}}{A character string giving the name(s) of the data.}
\item{\code{check.data}}{A \code{data.frame} containing the residuals for plotting.}
}
\description{
Projected covariance measure test for conditional mean independence
}
\details{
The projected covariance measure test tests whether the conditional
mean of Y given X and Z depends on X.
}
\examples{
n <- 150
X <- matrix(rnorm(2 * n), ncol = 2)
colnames(X) <- c("X1", "X2")
Z <- matrix(rnorm(2 * n), ncol = 2)
colnames(Z) <- c("Z1", "Z2")
Y <- rnorm(n) # X[, 2] + Z[, 2] + rnorm(1e3)
(pcm1 <- pcm(Y, X, Z))

}
\references{
Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2022). The
Projected Covariance Measure for assumption-lean variable significance
testing. arXiv preprint. \doi{10.48550/arXiv.2211.02039}
}
