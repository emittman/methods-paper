#Bayesian FDR

str(mcmc1)

#expected FDR approx. sum_g Pr(H_0g|Y)

library(plyr)

FDR <- ldply(1:10, function(j){
  sim <- readRDS(paste("sim_data/sim_", j, sep=""))
  true <- sim$truth$beta
  products <- lapply(C, function(hyp) hyp%*%t(true)>0)
  eval_hyp <- lapply(products, function(product) apply(product, 2, min))
  
  mcmc <- readRDS(paste("saved_mcmc/mcmc_sim_",j, sep=""))
  
  fdr <- ldply(1:6, function(i){
    hyp <- mcmc$summaries$probs[i,]
    ord  <- order(hyp, decreasing=TRUE)
    efdr <- cumsum(1-hyp[ord])/(1:length(hyp))
    fdr  <- cumsum(1 - eval_hyp[[i]][ord])/(1:length(hyp))
    data.frame(sim=j, iter = 1:length(hyp), hypothesis = names(eval_hyp)[i], efdr, fdr)
  })
})

filter(FDR, hypothesis == "high_mean") %>%
  ggplot(aes(efdr, fdr)) + geom_line() + facet_wrap(~sim) + geom_abline(slope=1, lty=2)

#very slow
ddply(FDR, .(iter, hypothesis), summarise,
      efdr = mean(efdr),
      fdr = mean(fdr)) %>%
filter(efdr<.2) %>%
ggplot(aes(efdr, fdr)) + 
  geom_line() + 
  geom_abline(slope=1, lty=2) + 
  facet_wrap(~hypothesis, scales="free")
