<!-- badges: start -->
[![R-CMD-check](https://github.com/LucasKook/comet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/comet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Covariance measure tests

The Generalised [1] and Projected [2] Covariance Measure tests (GCM, PCM) can be
used to test conditional independence between a real-valued response $Y$ and
features/modalities $X$ given additional features/modalities $Z$ in an
algorithm-agnostic way. The `comets` R package implements these covariance
measure tests (COMETs) with a user-friendly interface which allows the user to
use any sufficiently predictive supervised learning algorithm of their choosing.
The default is to use random forests implemented in `ranger` for all regressions.

Here, we showcase how to use `comets` with a simple example. More elaborate
examples including conditional variable significance testing and modality
selection on real-world data can be found in [3].

```r
set.seed(1)
n <- 300
X <- matrix(rnorm(2 * n), ncol = 2)
colnames(X) <- c("X1", "X2")
Z <- matrix(rnorm(2 * n), ncol = 2)
colnames(Z) <- c("Z1", "Z2")
Y <- X[, 1] + Z[, 2] + rnorm(n)
(GCM <- gcm(Y, X, Z)) # plot(GCM)
```

The output for the GCM test is shown below.

```
#	  Generalized covariance measure test
#
# data:  gcm(Y = Y, X = X, Z = Z)
# X-squared = 87.974, df = 2, p-value < 2.2e-16
# alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0
```

The PCM test can be run likewise.

```r
(PCM <- pcm(Y, X, Z)) # plot(PCM)
```

The outut is shown below.

```
#   Projected covariance measure test
#
# data:  pcm(Y = Y, X = X, Z = Z)
# Z = 5.5807, p-value = 1.198e-08
# alternative hypothesis: true E[Y | X, Z] is not equal to E[Y | Z]
```

The `comets` package contains an alternative formula-based interface,
in which $H_0 : Y \perp\hspace{-5pt}\perp X \mid Z$ can be supplied
as `Y ~ X | Z` with a corresponding `data` argument. This interface
is implemented in `comet()` and shown below.

```r
dat <- data.frame(Y = Y, X, Z)
comet(Y ~ X1 + X2 | Z1 + Z2, data = dat, test = "gcm")
```

Different regression methods can supplied for both GCM and PCM tests
using the `reg_*` arguments (for instance, `reg_YonZ` in `gcm()` for
the regression of $Y$ on $Z$). Pre-implemented regressions are `"rf"`
for random forests and `"lasso"` for cross-validated $L_1$-penalized
regression. Custom regression functions can be supplied as character
strings or functions, require a `predict()` method and the following
structure:

```
my_regression <- function(y, x, ...) {
  ret <- <run the regression>
  class(ret) <- "my_regression"
  ret
}

predict.my_regression <- function(object, data, ...) {
  <run the prediction rountine>
}
```

The input `y` and `x` and `data` are vector and matrix-valued. The output of 
`predict.my_regression()` should be a vector of `NROW(data)`.

# Replication materials

All results in [3] can be reproduced by running `make all` in `./inst`
after downloading all required data sets available [here](TODO).

# References

[1] Rajen D. Shah, Jonas Peters "The hardness of conditional independence
testing and the generalised covariance measure," The Annals of Statistics,
48(3), 1514-1538. [doi:10.1214/19-aos1857](https://doi.org/10.1214/19-aos1857)

[2] Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2022). The
Projected Covariance Measure for assumption-lean variable significance testing.
arXiv preprint.
[doi:10.48550/arXiv.2211.02039](https://doi.org/10.48550/arXiv.2211.02039)

[3] Kook, L. & Lundborg A. R. (2024). Algorithm-agnostic significance testing in
supervised learning with multimodal data. 
