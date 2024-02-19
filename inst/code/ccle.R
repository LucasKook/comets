### Pre-process CCLE data
### LK 2024

bpath <- "inst/data/ccle"
set.seed(2410)
do_screen <- TRUE

# DEPs --------------------------------------------------------------------

library("tidyverse")
library("survival")
devtools::load_all()

# Read --------------------------------------------------------------------

cmp <- "PLX4720"
nrep <- 10

# expression <- read_table("inst/data/ccle/expression.txt")
mutations <- read_table("inst/data/ccle/mutation.txt") |>
  filter(Name %in% grep("_MUT$", Name, value = TRUE))
response <- read_csv("inst/data/ccle/response.csv") |>
  filter(Compound == cmp)

idx <- match(response$`CCLE Cell Line Name`, colnames(mutations))

resp <- response |> filter(!is.na(idx))
mut <- mutations[, c(1, idx[!is.na(idx)])] |>
  column_to_rownames("Name")

cors <- apply(mut, 1, \(muts) {
  tc <- cor(resp$Amax, muts)
  ifelse(!is.na(tc), tc, 0)
})

screen <- if (do_screen) cors[abs(cors) > 0.05 & !is.na(cors)] else
  cors[!is.na(cors)]
dat <- data.frame(Y = resp$Amax, t(mut[names(screen),]))

# Run ---------------------------------------------------------------------

response <- "Y"
tgt <- paste0(c("BRAF.V600E", "BRAF.MC", "HIP1", "FLT3", "CDC42BPA",
                "THBS3", "DNMT1", "PRKD1", "PIP5K1A", "MAP3K5"), "_MUT")
muts <- setdiff(colnames(dat), c(response, tgt))

pb <- txtProgressBar(min = 0, max = length(tgt), style = 3)
out <- lapply(seq_along(tgt), \(tmut) {
  setTxtProgressBar(pb, tmut)
  pcm <- pcm(dat[, response], dat[, tgt[tmut]], dat[, c(tgt[-tmut], muts)],
             rep = nrep)
  gcm <- gcm(dat[, response], dat[, tgt[tmut]], dat[, c(tgt[-tmut], muts)])
  naive <- cor.test(x = dat[, tgt[tmut]], y = dat[, response])
  data.frame(gene = tgt[tmut],
             pcm = pcm$p.value,
             gcm = gcm$p.value,
             naive = naive$p.value)
})

### Tabulate results
(res <- bind_rows(out))
knitr::kable(t(res[, c("gcm", "pcm")]), col.names = res$gene,
             format = "latex", booktabs = TRUE, digits = 3)

### Save
write.csv(res, "inst/results/ccle-results.csv")
