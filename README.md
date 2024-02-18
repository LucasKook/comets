<!-- badges: start -->
[![R-CMD-check](https://github.com/LucasKook/comet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/comet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Covariance measure tests

The Generalised [1] and Projected [2] Covariance Measure tests (GCM, PCM) can be
used to test conditional independence between a real-valued response $Y$ and
features/modalities $X$ given additional features/modalities $Z$ in an
algorithm-agnostic way. The `comet` R package implements these covariance
measure tests (COMETs) with a user-friendly interface which allows the user to
use any sufficiently predictive supervised learning algorithm of their choosing.

Here, we showcase how to use `comet` with a simple example. More elaborate
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

```
#	  Generalized covariance measure test
#
# data:  gcm(Y = Y, X = X, Z = Z)
# X-squared = 87.974, df = 2, p-value < 2.2e-16
# alternative hypothesis: true E[cov(Y, X | Z)] is not equal to 0
```

```r
(PCM <- pcm(Y, X, Z)) # plot(PCM)
```

```
#   Projected covariance measure test
#
# data:  pcm(Y = Y, X = X, Z = Z)
# Z = 5.5807, p-value = 1.198e-08
# alternative hypothesis: true E[Y | X, Z] is not equal to E[Y | Z]
```

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
