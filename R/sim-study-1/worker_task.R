# i = identify the scenario 
# code below produce data, analysis and summary for ONE scenario

ss  = Sys.getenv("SLURM_ARRAY_TASK_ID") 
i = as.numeric(ss) + 1

# 0 ) libraries and functions
pkgs <- c('dplyr', 'cudarpackage', 'drake')
lapply(pkgs, library,  character.only = TRUE, quietly=TRUE )

# inputs

load('data/my_analyses.RData')
source('functions.R')
sim_name <- paste("sim_",i,sep="")
loadd(sim_name)
loadd(X)
this_mcmc <- run_mcmc(get(sim_name), X)
saveRDS(this_mcmc, file=paste("mcmc_sim_",i,sep=""))


