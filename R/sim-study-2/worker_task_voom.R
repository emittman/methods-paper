# i = identify the scenario 
# code below produce data, analysis and summary for ONE scenario

ss  = Sys.getenv("SLURM_ARRAY_TASK_ID") 
i = as.numeric(ss) + 1

# 0 ) libraries and functions
pkgs <- c('dplyr', 'cudarpackage', 'drake', 'edgeR', 'limma')
lapply(pkgs, library,  character.only = TRUE, quietly=TRUE )

# inputs

load('data/my_analyses.RData')
source('functions.R')
sim_name <- paste("sim_",i,sep="")
loadd()
this_mcmc <- run_mcmc(get(sim_name), X, voom=TRUE)
saveRDS(this_mcmc, file=paste("saved_mcmc/mcmc_voom_sim_",i,sep=""))