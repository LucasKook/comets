<!-- badges: start -->
[![R-CMD-check](https://github.com/LucasKook/comets/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/comets/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Covariance Measure Tests (COMETs) in R <img src='inst/comets-pkg.png' align="right" height="138.5" />

The Generalised [1], Projected [2], weighted generalised [3], kernel
generalised [4] Covariance Measure tests (GCM, PCM, wGCM, kGCM tests) can be
used to test conditional independence between a real-valued response $Y$ and
features/modalities $X$ given additional features/modalities $Z$ using any
sufficiently predictive supervised learning algorithms. An extension of the GCM
to censored and ordinal responses was proposed in [5] and is implemented with
survival regression methods and probability forests. The `comets` R package
implements these covariance measure tests (COMETs) with a user-friendly
interface which allows the user to use any sufficiently predictive supervised
learning algorithm of their choosing. The default is to use random forests
implemented in `ranger` for all regressions. A Python version of this package is
available
[here](https://github.com/shimenghuang/pycomets).

Here, we showcase how to use `comets` with a simple example in which $Y$ is not
independent of $X$ given $Z$. More elaborate examples including conditional
variable significance testing and modality selection on real-world data can be
found in [6].


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

The output for the GCM test, which fails to reject the null hypothesis of
conditional independence in this example, is shown below. The residuals for the
$Y$ on $Z$ and $X$ on $Z$ regressions can be investigated by calling `plot(GCM)`
(not shown here).


```
## 
## 	Generalised covariance measure test
## 
## data:  gcm(Y = Y, X = X, Z = Z)
## X-squared = 3.1556, df = 2, p-value = 0.2064
## alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0
```

The PCM test can be run likewise.


``` r
PCM <- pcm(Y, X, Z) # plot(PCM)
```

The output is shown below: The PCM test correctly rejects the null hypothesis of
conditional independence in this example.


```
## 
## 	Projected covariance measure test
## 
## data:  pcm(Y = Y, X = X, Z = Z)
## Z = 4.8727, p-value = 5.505e-07
## alternative hypothesis: true E[Y | X, Z] is not equal to E[Y | Z]
```

The `comets` package contains an alternative formula-based interface, in which
$H_0 : Y \perp\hspace{-5pt}\perp X \mid Z$ can be supplied as `Y ~ X | Z` with a
corresponding `data` argument. This interface is implemented in `comets()` and
shown below.


``` r
dat <- data.frame(Y = Y, X, Z)
comets(Y ~ X1 + X2 | Z1 + Z2, data = dat, test = "gcm")
```

```
## 
## 	Generalised covariance measure test
## 
## data:  comets(formula = Y ~ X1 + X2 | Z1 + Z2, data = dat, test = "gcm")
## X-squared = 3.2816, df = 2, p-value = 0.1938
## alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0
```

# Specifying regression methods

Different regression methods can supplied for both GCM and PCM tests using the
`reg_*` arguments (for instance, `reg_YonZ` in `gcm()` for the regression of $Y$
on $Z$). Pre-implemented regressions are `"rf"` for random forests and `"lasso"`
for cross-validated $L_1$-penalized regression. Custom regression functions can
be supplied as character strings or functions, require a `residual()` (GCM and
PCM) or `predict()` (PCM only) method and the following structure:

```
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
```

The input `y` and `x` and `data` are vector and matrix-valued. The output of
`predict.my_regression()` should be a vector of length `NROW(data)`.

# Usage example: Survival response

For survival responses, `comets` offers the TRAM-GCM test [5] and supports
parametric and semiparametric survival models from the `survival` package, as
well as random survival forests from `ranger`. As an example, we test whether
survival is independent of sex given age in the `cancer` dataset, once using a
Cox model and once using a random survival forest. Both tests agree to reject
the null hypothesis at conventional significance levels.


``` r
library("survival")
data("cancer", package = "survival")
cancer$surv <- with(cancer, Surv(time, status == 2))
comets(surv ~ sex | age, data = cancer, reg_YonZ = "cox")
```

```
## 
## 	TRAM-GCM test
## 
## data:  comets(formula = surv ~ sex | age, data = cancer, reg_YonZ = "cox")
## X-squared = 7.6583, df = 1, p-value = 0.005651
## alternative hypothesis: true E[cov(R, X | Z)] is not equal to 0
```

``` r
comets(surv ~ sex | age, data = cancer, reg_YonZ = "survforest")
```

```
## 
## 	TRAM-GCM test
## 
## data:  comets(formula = surv ~ sex | age, data = cancer, reg_YonZ = "survforest")
## X-squared = 14.824, df = 1, p-value = 0.000118
## alternative hypothesis: true E[cov(R, X | Z)] is not equal to 0
```

# Usage example: Multivariate response

The GCM test also supports multivariate responses. Continuing the example from
above, we generate a bivariate response $Y$. Internally, since $Y$ and $X$ are
both two dimensional, four random forest regressions are performed. Advanced
usage with the `multivariate` argument also allows the specification of
multivariate regression models (this option is experimental). 


``` r
bivY <- cbind(Y, 0.5 * X[, 1] + Z[, 1] + rnorm(n))
gcm(bivY, X, Z)
```

```
## 
## 	Generalised covariance measure test
## 
## data:  gcm(Y = bivY, X = X, Z = Z)
## X-squared = 40.301, df = 4, p-value = 3.75e-08
## alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0
```

# Usage example: Ordinal response

For ordered factor responses, the random forest regressions (`"rf"` and
`"tuned_rf"`) estimate the conditional distribution of $Y$ given $Z$ via a
probability forest. The argument `ordinal_residuals` controls which residuals
are computed from it: `"probability"` (default) gives the $K - 1$ dimensional
residuals $1(Y = k) - P(Y = k \mid Z)$, which have power against general
alternatives, `"mean"` gives $Y - E[Y \mid Z]$ using integer scores, and
`"score"` gives the score residual $P(Y' < Y \mid Z) - P(Y' > Y \mid Z)$ of a
proportional odds model whose transformation function is estimated by the
forest. With `"score"`, the GCM test is the TRAM-GCM test [5] for ordinal
responses. Below, we simulate an ordinal response with four levels that
depends on $X$ given $Z$ via a proportional odds model.


``` r
set.seed(1)
n <- 300
Z <- rnorm(n)
X <- Z + rnorm(n)
Y <- cut(Z + 0.3 * X + rlogis(n),
  breaks = c(-Inf, -1, 0, 1, Inf),
  labels = c("none", "mild", "moderate", "severe"), ordered_result = TRUE
)
dord <- data.frame(Y, X, Z)
comets(Y ~ X | Z, data = dord)
```

```
## 
## 	Generalised covariance measure test
## 
## data:  comets(formula = Y ~ X | Z, data = dord)
## X-squared = 8.849, df = 3, p-value = 0.03137
## alternative hypothesis: true E[cov(1(Y = k), X | Z)] is not equal to 0
```

``` r
comets(Y ~ X | Z, data = dord, args_YonZ = list(ordinal_residuals = "score"))
```

```
## 
## 	TRAM-GCM test
## 
## data:  comets(formula = Y ~ X | Z, data = dord, args_YonZ = list(ordinal_residuals = "score"))
## X-squared = 10.161, df = 1, p-value = 0.001434
## alternative hypothesis: true E[cov(R, X | Z)] is not equal to 0
```

Both tests reject the null hypothesis at the 5% level. The TRAM-GCM test based
on the one-dimensional score residual yields a smaller p-value, since the data
were simulated from a proportional odds model.

# Installation

The development version of `comets` can be installed using:

```r
# install.packages("remotes")
remotes::install_github("LucasKook/comets")
```

A stable version of `comets` can be installed from CRAN via:

```r
install.packages("comets")
```

# Replication materials

All results in [3] can be reproduced by running `make all` in `./inst` after
downloading all required data from the 
[zenodo repository](https://zenodo.org/doi/10.5281/zenodo.10689553).
The scripts for reproducing the results manually can be found in `./inst/code/`
for the CCLE data (`ccle.R`), TCGA data (`multiomics.R`) and MIMIC data
(`mimic.R`).

# Citation

Please cite the `comets` package as

```
@article{10.1093/bib/bbae475,
title={{Algorithm-Agnostic Significance Testing in Supervised Learning With Multimodal Data}},
author={Lucas Kook and Anton Rask Lundborg},
year={2024},
journal={Briefings in Bioinformatics},
volume={25},
number={6},
doi={10.1093/bib/bbae475},
}
```

# References

[1] Rajen D. Shah, Jonas Peters "The hardness of conditional independence
testing and the generalised covariance measure," The Annals of Statistics,
48(3), 1514-1538. [doi:10.1214/19-aos1857](https://doi.org/10.1214/19-aos1857)

[2] Lundborg, A. R., Kim, I., Shah, R. D., & Samworth, R. J. (2024). The
Projected Covariance Measure for assumption-lean variable significance testing.
The Annals of Statistics, 52(6), 2851-2878.
[doi:10.1214/24-AOS2447](https://doi.org/10.1214/24-AOS2447)

[3] Scheidegger, C., Hörrmann, J., & Bühlmann, P. (2022). The weighted
generalised covariance measure. Journal of Machine Learning Research, 23(273),
1-68.

[4] Fernández, T., & Rivera, N. (2024). A general framework for the analysis of
kernel-based tests. Journal of Machine Learning Research, 25(95), 1-40.

[5] Kook, L., Saengkyongam, S., Lundborg, A. R., Hothorn, T., & Peters, J.
(2025). Model-based causal feature selection for general response types.
Journal of the American Statistical Association, 120(550), 1090-1101.
[doi:10.1080/01621459.2024.2395588](https://doi.org/10.1080/01621459.2024.2395588)

[6] Kook, L. & Lundborg A. R. (2024). Algorithm-agnostic significance
testing in supervised learning with multimodal data. 
Briefings in Bioinformatics 25(6) 2024. 
[doi:10.1093/bib/bbae475](https://doi.org/10.1093/bib/bbae475)

