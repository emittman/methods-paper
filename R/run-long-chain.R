#long run, 1000000 iterations

source("workflow_fun.R")

seed <- 138762371
setup <- initialize_chain(seed, "stickBreaking")

saveRDS(setup, file="init-for-long-chain.rds")

s <- sample_bnp_model(setup)

saveRDS(s, file="long-chain-results.rds")