### Pre-process multiomics data
### LK 2024

bpath <- "inst/data/multiomics"
set.seed(2410)

if (!dir.exists("inst/results"))
  dir.create("inst/results")

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
n_mir <- grep("^mir", colnames(dat), value = TRUE)
n_rna <- grep("^rna", colnames(dat), value = TRUE)
n_meth <- grep("^meth", colnames(dat), value = TRUE)

nrep <- 10

(t1 <- pcm(dat[, response][, 2], dat[, n_rna], dat[, c(n_mir, n_meth)],
           rep = nrep))
(t2 <- pcm(dat[, response][, 2], dat[, n_mir], dat[, c(n_rna, n_meth)],
           rep = nrep))
(t3 <- pcm(dat[, response][, 2], dat[, n_meth], dat[, c(n_rna, n_mir)],
           rep = nrep))

(t4 <- pcm(dat[, response][, 2], dat[, n_rna], dat[, c(n_mir, n_meth)],
           rep = nrep, reg = "lasso"))
(t5 <- pcm(dat[, response][, 2], dat[, n_mir], dat[, c(n_rna, n_meth)],
           rep = nrep, reg = "lasso"))
(t6 <- pcm(dat[, response][, 2], dat[, n_meth], dat[, c(n_rna, n_mir)],
           rep = nrep, reg = "lasso"))
