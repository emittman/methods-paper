source("functions.R")
library(cudarpackage)
load("data/counts_1.RData")
load("data/design.RData")
load("data/truth.RData")


s <- run_mcmc(data.matrix(y), X)

saveRDS(s, file="mcmc_1.rds")