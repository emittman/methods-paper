#long run, 300000 iterations
library(cudarpackage)
source("workflow_fun.R")

seed <- 138762371
setup <- initialize_chain(seed, "stickBreaking", 10000, 10000, 300000)

saveRDS(setup, file="init-for-long-chain.rds")

s <- sample_bnp_model(setup)

saveRDS(s, file="long-chain-results.rds")