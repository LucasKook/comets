### Timings PCM/GCM with random forest
### LK 2024

# Dependencies ------------------------------------------------------------

library("comets")
library("tidyverse")
save <- TRUE

# Functions ---------------------------------------------------------------

dgp <- function(n = 50, p = 10) {
  Z <- matrix(rlogis(n * p), nrow = n, ncol = p)
  X <- Z + matrix(rlogis(n * p), nrow = n, ncol = p)
  Y <- rowSums(Z) + rlogis(n)
  list(Y = Y, X = X, Z = Z)
}

# Parameters --------------------------------------------------------------

niter <- 10
ns <- 50 * 2^(0:8)
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
  geom_boxplot(position = position_dodge(0), outlier.size = rel(0.1), size = rel(0.5)) +
  facet_wrap(~ p, labeller = label_bquote(d==.(p)), nrow = 2) +
  labs(x = "sample size", y = "elapsed time (in seconds)", color = "COMET") +
  theme_bw() +
  scale_y_log10() +
  scale_x_log10() +
  theme(text = element_text(size = 13.5), legend.position = "top") +
  scale_color_brewer(palette = "Dark2")

# Save --------------------------------------------------------------------

if (save) {
  if (!dir.exists("inst/results"))
    dir.create("inst/results", recursive = TRUE)

  write_csv(res, "inst/results/timings.csv")
  ggsave("inst/figures/timings.pdf", height = 5, width = 7.5, scale = 0.8)
}
