### Pre-process multiomics data
### LK 2024

bpath <- "inst/data/multiomics"
set.seed(2410)

# DEPs --------------------------------------------------------------------

library("tidyverse")
library("survival")
devtools::load_all()

# Read --------------------------------------------------------------------

rna <- read_tsv(file.path(bpath, "rna.tsv"))
meth <- read_tsv(file.path(bpath, "meth.tsv"))
mir <- read_tsv(file.path(bpath, "mir.tsv"))
surv <- read_tsv(file.path(bpath, "survival.tsv")) %>%
  mutate(surv = Surv(days, event))

rna_idx <- match(rna$Samples, surv$Samples)
mir_idx <- match(mir$Samples, surv$Samples)
meth_idx <- match(meth$Samples, surv$Samples)

stopifnot(all(rna_idx == mir_idx))
stopifnot(all(meth_idx == mir_idx))

dat <- data.frame(
  surv = surv$surv[rna_idx],
  rna = scale(as.matrix(rna[, -1])),
  mir = scale(as.matrix(mir[, -1])),
  meth = scale(as.matrix(meth[, -1]))
)

rm(list = c("mir", "rna", "surv", "meth"))

# Run ---------------------------------------------------------------------

response <- "surv"
names_mir <- grep("^mir", colnames(dat), value = TRUE)
names_rna <- grep("^rna", colnames(dat), value = TRUE)
names_meth <- grep("^meth", colnames(dat), value = TRUE)

rargs <- list(mtry = NULL, max.depth = NULL, probability = TRUE)
tmtry <- NULL
nrep <- 10

(t1 <- pcm(dat[, response][, 2], dat[, names_rna], dat[, c(names_mir, names_meth)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))
(t2 <- pcm(dat[, response][, 2], dat[, names_mir], dat[, c(names_rna, names_meth)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))
(t3 <- pcm(dat[, response][, 2], dat[, names_meth], dat[, c(names_rna, names_mir)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))

(t4 <- pcm(dat[, response][, 2], dat[, names_rna], dat[, c(names_mir, names_meth)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, reg = "pcm_lasso"))
(t5 <- pcm(dat[, response][, 2], dat[, names_mir], dat[, c(names_rna, names_meth)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, reg = "pcm_lasso"))
(t6 <- pcm(dat[, response][, 2], dat[, names_meth], dat[, c(names_rna, names_mir)],
           mtry = tmtry, rep = nrep, est_vhat = TRUE, reg = "pcm_lasso"))

### Don't run
# (g1 <- ranger_gcm(dat[, response][, 2], dat[, names_rna], dat[, c(names_mir, names_meth)],
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))
# (g2 <- ranger_gcm(dat[, response][, 2], dat[, names_mir], dat[, c(names_rna, names_meth)],
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))
# (g3 <- ranger_gcm(dat[, response][, 2], dat[, names_meth], dat[, c(names_rna, names_mir)],
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))

# screen_rna <- names(which(sapply(names_rna, \(tx) cor(dat[, response][, 2], dat[, tx])) > 0.2))
# screen_mir <- names(which(sapply(names_mir, \(tx) cor(dat[, response][, 2], dat[, tx])) > 0.2))
# screen_meth <- names(which(sapply(names_meth, \(tx) cor(dat[, response][, 2], dat[, tx])) > 0.2))
#
# (t1 <- pcm(dat[, response][, 2], dat[, screen_rna], dat[, c(screen_mir, screen_meth)],
#            mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))
# (t2 <- pcm(dat[, response][, 2], dat[, screen_mir], dat[, c(screen_rna, screen_meth)],
#            mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))
# (t3 <- pcm(dat[, response][, 2], dat[, screen_meth], dat[, c(screen_rna, screen_mir)],
#            mtry = tmtry, rep = nrep, est_vhat = TRUE, ghat_args = rargs))
# (g1 <- ranger_gcm(dat[, response][, 2], as.matrix(dat[, screen_rna]),
#                   as.matrix(dat[, c(screen_mir, screen_meth)]),
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))
# (g2 <- ranger_gcm(dat[, response][, 2], as.matrix(dat[, screen_mir]),
#                   as.matrix(dat[, c(screen_rna, screen_meth)]),
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))
# (g3 <- ranger_gcm(dat[, response][, 2], as.matrix(dat[, screen_meth]),
#                   as.matrix(dat[, c(screen_rna, screen_mir)]),
#                   mtry = rargs$mtry, max.depth = rargs$max.depth))

# SurvGCM -----------------------------------------------------------------

# dmeth <- dat[, c(response, screen_meth)]
# pb <- txtProgressBar(0, length(screen_meth), style = 3)
# out <- lapply(seq_along(screen_meth), \(tmeth) {
#   setTxtProgressBar(pb, tmeth)
#   fm <- reformulate(screen_meth[-tmeth], response)
#   mY <- survforest(fm, dmeth)
#   tst <- gcm(mY, reformulate(screen_meth[-tmeth], screen_meth[tmeth]), data = dmeth)
#   mY2 <- .ranger(reformulate(screen_meth[-tmeth], "surv[, 2]"), dmeth)
#   tst2 <- gcm(mY2, reformulate(screen_meth[-tmeth], screen_meth[tmeth]), data = dmeth)
#   # mY2 <- coxph(fm, dmeth)
#   # tst2 <- gcm(mY2, reformulate(screen_meth[-tmeth], screen_meth[tmeth]), data = dmeth)
#   data.frame(meth = screen_meth[tmeth], survgcm = tst$p.value, gcm = tst2$p.value)
# })
# bind_rows(out)
#
# ggplot(bind_rows(out), aes(x = gcm, y = coxph, color = gcm < 0.05 & coxph < 0.05)) +
#   geom_point(show.legend = FALSE) +
#   theme_bw() +
#   scale_x_log10() +
#   scale_y_log10() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(x = "Survival GCM", y = "Binary GCM")
