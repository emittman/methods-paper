#mse


load("../data/heterosis_long.RData")
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

library(limma)
library(plyr)

mse <- ldply(1:10, function(i){
  sim <- readRDS(paste("sim_data/sim_",i,sep=""))
  y <- sim$y
  
  limma_fit <- lmFit(y, design = X)
  limma_fit <- eBayes(limma_fit)
  
  beta_limma <- coefficients(limma_fit)
  beta_true <- sim$truth$beta
  beta_bnp <- t(readRDS(paste("saved_mcmc/mcmc_sim_",i,sep=""))$summaries$means_betas)
  
  se_limma <- (beta_limma - beta_true)^2
  mse_limma <- data.frame(t(apply(se_limma, 2, mean)))
  se_bnp <- (beta_bnp - beta_true)^2
  mse_bnp <- data.frame(t(apply(se_bnp, 2, mean)^2))
  names(mse_limma) <- sapply(1:5, function(p) paste("beta[",p,"]",sep=""))
  names(mse_bnp) <- sapply(1:5, function(p) paste("beta[",p,"]",sep=""))
  out <- rbind(mse_limma, mse_bnp)
  out$type <- c("limma","bnp")
  out$sim <- i
  out
})

save(mse, file="summaries/mse.RData")
load("summaries/mse.RData")
library(ggplot2)
library(tidyr)
gather(mse, key=par, value=value, -type, -sim) %>%
  ggplot(aes(factor(par), value, color=type)) +
  geom_boxplot(position=position_dodge()) +
  theme_bw()
