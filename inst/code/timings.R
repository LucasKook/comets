### Timings PCM/GCM with random forest
### LK 2024

# Dependencies ------------------------------------------------------------

library("comets")
library("tidyverse")

# Functions ---------------------------------------------------------------

dgp <- function(n = 50, p = 10) {
  Z <- matrix(rlogis(n * p), nrow = n, ncol = p)
  X <- Z + matrix(rlogis(n * p), nrow = n, ncol = p)
  Y <- rowSums(Z) + rlogis(n)
  list(Y = Y, X = X, Z = Z)
}

# Parameters --------------------------------------------------------------

niter <- 20
ns <- 50 * 2^(0:6)
ps <- 1 * 2^(0:5)

# Run ---------------------------------------------------------------------

res <- lapply(ns, \(tn) {
  cat("Running n =", tn, "\n")
  lapply(ps, \(tp) {
    cat("Running p =", tp, "\n")
    lapply(1:niter, \(iter) {
      d <- dgp(tn, tp)
      gcm <- system.time(gcm(d$Y, d$X, d$Z, coin = TRUE))["elapsed"]
      pcm <- system.time(pcm(d$Y, d$X, d$Z, coin = TRUE))["elapsed"]
      data.frame(iter = iter, n = tn, p = tp, GCM = gcm, PCM = pcm)
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows()

# Visualize ---------------------------------------------------------------

pdat <- res |>
  pivot_longer(cols = c("GCM", "PCM"), names_to = "test", values_to = "time")

ggplot(pdat, aes(x = n, y = time, color = test, group = interaction(n, test))) +
  geom_boxplot() +
  facet_wrap(~ p, labeller = label_bquote(p==.(p))) +
  labs(x = "sample size", y = "elapsed time (in seconds)", color = "COMET") +
  theme_bw() +
  scale_y_log10() +
  scale_x_log10()
