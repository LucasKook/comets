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
