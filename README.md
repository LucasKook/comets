<!-- badges: start -->

[![R-CMD-check](https://github.com/LucasKook/comets/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/comets/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Covariance Measure Tests (COMETs) in R <img src='inst/comets-pkg.png' align="right" height="138.5" />

The Generalised \[1\] and Projected \[2\] Covariance Measure tests (GCM,
PCM) can be used to test conditional independence between a real-valued
response *Y* and features/modalities *X* given additional
features/modalities *Z* using any sufficiently predictive supervised
learning algorithms. The `comets` R package implements these covariance
measure tests (COMETs) with a user-friendly interface which allows the
user to use any sufficiently predictive supervised learning algorithm of
their choosing. The default is to use random forests implemented in
`ranger` for all regressions. A Python version of this package is
available [here](https://github.com/shimenghuang/pycomets).

Here, we showcase how to use `comets` with a simple example in which *Y*
is not independent of *X* given *Z*. More elaborate examples including
conditional variable significance testing and modality selection on
real-world data can be found in \[3\].

``` r
set.seed(1)
n <- 300
X <- matrix(rnorm(2 * n), ncol = 2)
colnames(X) <- c("X1", "X2")
Z <- matrix(rnorm(2 * n), ncol = 2)
colnames(Z) <- c("Z1", "Z2")
Y <- X[, 1]^2 + Z[, 2] + rnorm(n)
GCM <- gcm(Y, X, Z) # plot(GCM)
```

The output for the GCM test, which fails to reject the null hypothesis
of conditional independence in this example, is shown below. The
residuals for the *Y* on *Z* and *X* on *Z* regressions can be
investigated by calling `plot(GCM)` (not shown here).

    ## 
    ##  Generalized covariance measure test
    ## 
    ## data:  gcm(Y = Y, X = X, Z = Z)
    ## X-squared = 2.8211, df = 2, p-value = 0.244
    ## alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0

The PCM test can be run likewise.

``` r
PCM <- pcm(Y, X, Z) # plot(PCM)
```

The output is shown below: The PCM test correctly rejects the null
hypothesis of conditional independence in this example.

    ## 
    ##  Projected covariance measure test
    ## 
    ## data:  pcm(Y = Y, X = X, Z = Z)
    ## Z = 4.7107, p-value = 1.235e-06
    ## alternative hypothesis: true E[Y | X, Z] is not equal to E[Y | Z]

The `comets` package contains an alternative formula-based interface, in
which $H_0 : Y \perp\hspace{-5pt}\perp X \mid Z$ can be supplied as
`Y ~ X | Z` with a corresponding `data` argument. This interface is
implemented in `comet()` and shown below.

``` r
dat <- data.frame(Y = Y, X, Z)
comet(Y ~ X1 + X2 | Z1 + Z2, data = dat, test = "gcm")
```

    ## 
    ##  Generalized covariance measure test
    ## 
    ## data:  test(Y = Y, X = X, Z = Z)
    ## X-squared = 3.2184, df = 2, p-value = 0.2
    ## alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0

Different regression methods can supplied for both GCM and PCM tests
using the `reg_*` arguments (for instance, `reg_YonZ` in `gcm()` for the
regression of *Y* on *Z*). Pre-implemented regressions are `"rf"` for
random forests and `"lasso"` for cross-validated
*L*<sub>1</sub>-penalized regression. Custom regression functions can be
supplied as character strings or functions, require a `residual()` (GCM
and PCM) or `predict()` (PCM only) method and the following structure:

    my_regression <- function(y, x, ...) {
      ret <- <run the regression>
      class(ret) <- "my_regression"
      ret
    }

    predict.my_regression <- function(object, data, ...) {
      <run the prediction routine>
    }

    residuals.my_regression <- function(object, response, data, ...) {
      <run the routine for computing residuals>
    }

The input `y` and `x` and `data` are vector and matrix-valued. The
output of `predict.my_regression()` should be a vector of length
`NROW(data)`.

# Installation

The development version of `comets` can be installed using:

``` r
# install.packages("remotes")
remotes::install_github("LucasKook/comets")
```

A stable version of `comets` can be installed from CRAN via:

``` r
install.packages("comets")
```

# Replication materials

All results in \[3\] can be reproduced by running `make all` in `./inst`
after downloading all required data from the [zenodo
repository](https://zenodo.org/doi/10.5281/zenodo.10689553). The scripts
for reproducing the results manually can be found in `./inst/code/` for
the CCLE data (`ccle.R`), TCGA data (`multiomics.R`) and MIMIC data
(`mimic.R`).

# References

\[1\] Rajen D. Shah, Jonas Peters “The hardness of conditional
independence testing and the generalised covariance measure,” The Annals
of Statistics, 48(3), 1514-1538.
[doi:10.1214/19-aos1857](https://doi.org/10.1214/19-aos1857)

\[2\] Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2022).
The Projected Covariance Measure for assumption-lean variable
significance testing. arXiv preprint.
[doi:10.48550/arXiv.2211.02039](https://doi.org/10.48550/arXiv.2211.02039)

\[3\] Kook, L. & Lundborg A. R. (2024). Algorithm-agnostic significance
testing in supervised learning with multimodal data. arXiv preprint.
[doi:10.48550/arXiv.2402.14416](https://doi.org/10.48550/arXiv.2402.14416)
