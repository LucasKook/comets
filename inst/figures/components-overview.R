### Figure 1
### LK 2024

set.seed(2410)

# DEPs --------------------------------------------------------------------

library("tidyverse")
library("mgcv")
devtools::load_all()
theme_set(theme_bw() + theme(text = element_text(size = 13.5)))

pd1 <- data.frame(id = seq_len(1e3), pvalue = -log10(p.adjust(runif(1e3)^2)))

ggplot(pd1, aes(x = id, y = pvalue, col = pvalue)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), color = "darkred", linetype = 2) +
  labs(x = element_blank(), y = parse(text = "-log[10]~p*'-'*value"),
       color = parse(text = "-log[10]~p*'-'*value")) +
  scale_color_viridis_c(begin = 0.1, end = 0.9, option = "C") +
  theme(legend.position = "top", axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("inst/figures/manhatten.pdf", height = 2, width = 3)

# GCM v PCM ---------------------------------------------------------------

f <- \(x) 1 + sin(3 * x^2)
g <- \(z) 1 + z^3
mu <- \(x, z) f(x) * g(z)

X <- runif(n <- 3e2, -1, 1)
Z <- runif(n, -1, 1)
Y <- mu(X, Z) + rnorm(n, sd = 0.05)
pd2 <- data.frame(Y = Y, X = X, Z = Z)

mY <- .ranger(Y ~ Z, data = pd2)
pd2$RY <- residuals.ranger(mY)
mX <- .ranger(X ~ Z, data = pd2)
pd2$RX <- residuals.ranger(mX)
pd2$TX <- (f(X) - mean(f(X))) * g(Z)
mTX <- .ranger(TX ~ Z, data = pd2)
pd2$RTX <- residuals.ranger(mTX)

ggplot(pd2 |> mutate(id = 1:n) |> pivot_longer(X:Z, names_to = "Predictor",
                                               values_to = "value"),
       aes(y = Y, x = value, color = Predictor, group = Predictor)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", se = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.box.spacing = unit(-0.2, "cm")) +
  labs(y = "Response Y", subtitle = "Y depends on X given Z",
       x = element_blank())

ggsave("inst/figures/pcmgcmdata.pdf", height = 3.3, width = 3.3)

ggplot(pd2, aes(x = X, y = RY)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, col = "gray80") +
  geom_smooth(method = "gam", se = FALSE, col = "#1E88E5", linetype = 2, linewidth = 1.8) +
  labs(x = "X", y = "Residual Y on Z",
       subtitle = "GCM (fails to reject)") +
  geom_label(aes(x = -0.5, y = 0.05,
                 label = paste0("R==", round(cor(RX, RY), 2))),
             inherit.aes = FALSE, parse = TRUE)

ggsave("inst/figures/gcm.pdf", height = 3.3, width = 3.3)

ggplot(pd2, aes(x = TX, y = RY)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, col = "gray80") +
  labs(x = "Transformation of X and Z", y = "Residual Y on Z",
       subtitle = "PCM (correctly rejects)") +
  geom_label(aes(x = 0.3, y = 0.1,
                 label = paste0("R==", round(cor(TX, RY), 2))),
             inherit.aes = FALSE, parse = TRUE) +
  xlim(-1.05, 1.05)

ggsave("inst/figures/pcm.pdf", height = 3.3, width = 3.3)

# Modality selection ------------------------------------------------------

pd3 <- data.frame(MSE1 = runif(3e1, min = 1, max = 3),
                  MSE2 = runif(3e1, min = 0.3, max = 2.5))
ggplot(pivot_longer(pd3, MSE1:MSE2), aes(x = name, y = value)) +
  geom_violin(width = 0.5) +
  geom_boxplot(width = 0.3) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3, width = 0.2) +
  labs(x = element_blank(), y = "MSPE") +
  scale_x_discrete(labels = c("MSE1" = "Y | X", "MSE2" = "Y | X, Z"))

ggsave("inst/figures/modsec.pdf", height = 2, width = 3)

# GCM ---------------------------------------------------------------------

X <- matrix(runif(3e2), ncol = 1)
colnames(X) <- c("X1")
Z <- matrix(runif(6e2), ncol = 2)
colnames(Z) <- c("Z1", "Z2")
Y <- exp(X[, 1]) + Z[, 2]^3 + exp(Z[, 1]) + rnorm(3e2)
(gcm1 <- gcm(Y, X, Z))
plot(gcm1) + labs(x = "Transformation of X and Z")

ggsave("inst/figures/simple-gcm.pdf", height = 2, width = 2.5)
