library(drake)
library(dplyr)
library(cudarpackage)

source("functions.R")

G <- 10000

#load Psamples
load("data/P-samples.RData")

#load heterosis data ... to get design matrix
load("../data/heterosis_long.RData")
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

my_data <- expand(plan(sim = generate_log_counts(G, X, Psamples)), 1:10)

my_analyses <- analyses(plan(mcmc = run_mcmc(..dataset.., X)
                          #,
                          #limma = fit_limma(..dataset..)
                          ), datasets = my_data)
#compute coverage of 95% credible intervals for beta_{gj}



make(my_data, jobs=3)

save(my_analyses, file="data/my_analyses.RData")
