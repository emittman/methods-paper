library(drake)

source("functions.R")
load("../data/heterosis_design.RData")

setup <- sim_study_2_setup()
G <- nrow(setup$betas)

sims <- expand(plan(sim = generate_counts(G, 10000, X, setup)), values = 1:10)

make(plan=sims)