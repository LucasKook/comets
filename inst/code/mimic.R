### MIMIC analysis GCM/PCM
### LK 2023

set.seed(2410)
nrep <- 5
red <- 111 # 98% variance

if (!dir.exists("inst/results"))
  dir.create("inst/results")

# Dependencies ------------------------------------------------------------

library("readr")
library("ranger")
library("ggplot2")
devtools::load_all()

# Read data ---------------------------------------------------------------

## Read meta data
bpath <- "inst/data/mimic"

if (file.exists(temb <- file.path(bpath, "emb.csv"))) {
  emb <- as.data.frame(read_csv(temb))
  emb$race <- as.factor(emb$race)
  emb$sex <- as.factor(emb$sex)
} else {
  library("reticulate")
  np <-import("numpy")
  meta <- read_csv(file.path(bpath, "mimic_cfm_train_meta.csv"))
  meta_test <- read_csv(file.path(bpath, "mimic_cfm_test_meta.csv"))
  meta$sex <- as.factor(meta$sex)
  meta_test$sex <- as.factor(meta_test$sex)
  ## race: Merge Asian, White, Black
  meta$race <- as.factor(meta$race)
  meta_test$race <- as.factor(meta_test$race)

  if(file.exists(file.path(bpath, "emb_svd.RDS"))) {

    emb <- np$load(file.path(bpath, "mimic_cfm_train_emb.npy"))
    emb_svd <- readRDS(file.path(bpath, "emb_svd.RDS"))

  } else {

    emb <- np$load(file.path(bpath, "mimic_cfm_train_emb.npy"))

    ## SVD
    emb_svd <- fast.svd(emb)
    saveRDS(emb_svd, file = file.path(bpath, "emb_svd.RDS"))
    plot(cumsum(emb_svd$d^2/sum(emb_svd$d^2)), type="b")

  }

  # create reduced embedding from SVD
  emb <- emb %*% emb_svd$v[, 1:red]
  emb_test <- tcrossprod(
    np$load(file.path(bpath, "mimic_cfm_test_emb.npy")), emb_svd$v[1:red,])

  ## Combine
  emb <- cbind(as.data.frame(emb), meta)
  emb_test <- cbind(as.data.frame(emb_test), meta_test)
  rm(meta, meta_test); gc()

  ## Create response
  name_resp <- "Pleural Effusion"
  emb$resp <- emb[[name_resp]]

  emb <- emb |> dplyr::select(
    dplyr::matches(name_resp), sex, age, race, dplyr::starts_with("V"))

  write_csv(emb, temb)
}

# Test --------------------------------------------------------------------

### Variables
nY <- "Pleural Effusion"
nX <- "race"
nC <- c("sex", "age")
nZ <- paste0("V", 1:red)

### Test function
run_tests <- function(splits = 20, max_size = 1e4, verbose = FALSE) {
  pb <- txtProgressBar(min = 0, max = splits, style = 3)
  n <- NROW(emb)
  folds <- sample(rep(1:splits, ceiling(n/splits)), n)
  lapply(seq_len(splits), \(iter) {
    setTxtProgressBar(pb, iter)
    idx <- which(folds == iter)
    if (length(idx) > max_size)
      idx <- idx[1:max_size]
    nemb <- emb[idx, ]

    ### Test resp _||_ race | emb, age, sex ### GCM/PCM
    gcm1 <- gcm(nemb[, nY], .mm(nX, nemb), .mm(c(nZ, nC), nemb))
    if (verbose)
      cat("\nGCM1 done")
    pcm1 <- pcm(nemb[, nY], .mm(nX, nemb), .mm(c(nZ, nC), nemb), rep = nrep)
    if (verbose)
      cat("\nPCM1 done")
    ### Test resp _||_ emb | race, age, sex ### GCM/PCM
    gcm2 <- gcm(nemb[, nY], .mm(nZ, nemb), .mm(c(nX, nC), nemb))
    if (verbose)
      cat("\nGCM2 done")
    pcm2 <- pcm(nemb[, nY], .mm(nZ, nemb), .mm(c(nX, nC), nemb), rep = nrep)
    if (verbose)
      cat("\nPCM2 done")
    ### Return
    list(GCM1 = gcm1, PCM1 = pcm1, GCM2 = gcm2, PCM2 = pcm2)
  })
}

### RUN on full data
print(full <- run_tests(1, nrow(emb), verbose = TRUE))
pvals <- c(
  "GCM1" = -log10(full[[1]]$GCM1$p.value),
  "PCM1" = -log10(full[[1]]$PCM1$p.value),
  "GCM2" = -pchisq(unname(full[[1]]$GCM2$statistic), log.p = TRUE,
                   df = 111, lower.tail = FALSE) / log(10),
  "PCM2" = -pnorm(unname(full[[1]]$PCM2$statistic), log.p = TRUE,
                  lower.tail = FALSE) / log(10)
)
cat("Negative log10 p-values\n")
print(pvals)
saveRDS(full, "inst/results/mimic-full.rds")

### RUN with splits
ms <- 150 * 4^(0:2)
nsplit <- 75
out <- dplyr::bind_rows(lapply(ms, \(nmax) {
  res <- run_tests(splits = nsplit, max_size = nmax)
  dplyr::bind_rows(lapply(res, \(x) unlist(lapply(x, \(y) y$p.value)))) |>
    dplyr::mutate(n = nmax)
}))

out |> tidyr::pivot_longer(-n) |>
  dplyr::mutate(hypothesis = dplyr::case_when(
    name %in% c("GCM1", "PCM1") ~ "PE~'_||_'~race~'|'~x*'-'*ray*','*~sex*','*~age",
    name %in% c("GCM2", "PCM2") ~ "PE~'_||_'~x*'-'*ray~'|'~race*','*~sex*','*~age"
  )) |>
  ggplot(aes(x = ordered(n), color = stringr::str_remove(name, "[0-9]"), y = -log10(value))) +
  geom_violin(width = 0.5, position = position_dodge(width = 0.5)) +
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.5), outlier.shape = NA) +
  ggbeeswarm::geom_quasirandom(width = 0.1, dodge.width = 0.5, alpha = 0.3, size = 0.5) +
  facet_wrap(~ hypothesis, labeller = label_parsed, scales = "free") +
  labs(x = "Sample size", y = parse(text = "-log[10]~p*'-'*value"), color = "COMET") +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = "darkred", linetype = 3) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "top")

readr::write_csv(out, "inst/results/mimic-pvals.csv")
ggsave("inst/figures/mimic-pval.pdf", height = 3.5, width = 6, scale = 0.85)
