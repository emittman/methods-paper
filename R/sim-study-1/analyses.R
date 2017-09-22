library(drake)
library(dplyr)
library(cudarpackage)

source("functions.R")

G <- 1000

mcmc <- readRDS("../arrayOutput/chain1_sb.rds")
Psamples <- mcmc$samples$P

load("../data/heterosis_long.RData")
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

my_data <- expand(plan(sim = generate_log_counts(G, X, Psamples)), 1:10)

my_analyses <- analyses(plan(mcmc = run_mcmc(..dataset.., X)
                          #,
                          #limma = fit_limma(..dataset..)
                          ), datasets = my_data)

make(my_data, jobs=1)
make(my_analyses, jobs=1)
