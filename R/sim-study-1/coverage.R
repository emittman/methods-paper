#coverage
library(plyr)

combined_coverage <- ldply(1:10, function(i){

  mcmcfile <- paste("saved_mcmc/mcmc_sim_",i,sep="")
  simfile <- paste("sim_data/sim_",i,sep="")
  mcmc <- readRDS(mcmcfile)
  sim <- readRDS(simfile)
  
  beta_samp <- mcmc$samples$beta
  beta_true <- sim$truth$beta
  
  # str(beta_samp)
  # str(beta_true)
  
  cis <- adply(beta_samp, c(1,3), function(beta_gp){
    data.frame(t(quantile(beta_gp, c(.05,.95))))
  })
  cis$truth <- as.numeric(t(beta_true))
  
  cis <- dplyr::mutate(cis, covered = truth < X95. & truth > X5.)
  
  coverage <- ddply(cis, .(v), summarise, coverage_rate = mean(covered))
  coverage$sim <- i
  coverage

})

saveRDS(combined_coverage, file="summaries/combined_coverage.rds")
combined_coverage <- readRDS("summaries/combined_coverage.rds")
ggplot(combined_coverage, aes(v, coverage_rate)) + geom_boxplot() + geom_hline(yintercept=.9, lty=2)

# mcmc1 <- readRDS("saved_mcmc/mcmc_sim_1")
# sim_1 <- readRDS("sim_data/sim_1")
# 
# beta_samp1 <- mcmc1$samples$beta
# beta_true1 <- sim_1$truth$beta
# 
# str(beta_samp1)
# str(beta_true1)
# 
# library(plyr)
# cis_1 <- adply(beta_samp1, c(1,3), function(beta_gp){
#   data.frame(t(quantile(beta_gp, c(.05,.95))))
# })
# cis_1$truth <- as.numeric(t(beta_true1))
# 
# cis_1 <- dplyr::mutate(cis_1, covered = truth < X95. & truth > X5.)
# 
# ddply(cis_1, .(v), summarise, coverage_rate = mean(covered))
# 
# id <- sample(1000, 100)
# dplyr::filter(cis_1, g %in% id) %>%
# ggplot(aes(x=g, ymin=X5., ymax=X95., y=truth, color=covered)) + geom_pointrange() + facet_wrap(~v, scales="free")

