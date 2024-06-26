% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{plm_equiv_test}
\alias{plm_equiv_test}
\title{Equivalence test for the parameter in a partially linear model}
\usage{
plm_equiv_test(Y, X, Z, from, to, ...)
}
\arguments{
\item{Y}{Vector or matrix of response values.}

\item{X}{Matrix or data.frame of covariates.}

\item{Z}{Matrix or data.frame of covariates.}

\item{from}{Lower bound of the equivalence margin}

\item{to}{Upper bound of the equivalence margin}

\item{...}{Further arguments passed to \code{\link{gcm}}}
}
\value{
Object of class '\code{gcm}' and '\code{htest}'
}
\description{
Equivalence test for the parameter in a partially linear model
}
\details{
The partially linear model postulates \deqn{Y = X \theta + g(Z) + \epsilon,}
and the target of inference is theta. The target is closely related to
the conditional covariance between Y and X given Z:
\deqn{\theta = E[cov(X, Y | Z)] / E[Var(X | Z)].} The equivalence test (based
on the GCM test) tests \eqn{H_0: \theta \not\in [{\tt from}, {\tt to}]} versus
\eqn{H_1: \theta \in [{\tt from}, {\tt to}]}. Y, X (and theta) can only be
one-dimensional. There are no restrictions on Z.
}
\examples{
n <- 150
X <- rnorm(n)
Z <- matrix(rnorm(2 * n), ncol = 2)
colnames(Z) <- c("Z1", "Z2")
Y <- X^2 + Z[, 2] + rnorm(n)
plm_equiv_test(Y, X, Z, from = -1, to = 1)
}
