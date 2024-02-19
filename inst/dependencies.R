# Install all dependencies for reproducing the results in the paper
# LK 2024

.libPaths(c("~/tutu/lib", .libPaths()))

### Install remotes for installing all other packages
install.packages("remotes", repos = "https://cloud.r-project.org")

### Install all other packages from CRAN
pkgs <- c("tidyverse", "ranger", "mlt", "sandwich", "glmnet", "tram",
          "reticulate", "mgcv", "survival")
remotes::install_cran(pkgs)

### Install comets
remotes::install_local(force = TRUE)
