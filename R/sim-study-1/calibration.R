#calibrated posterior probabilities?
#to do for each (sim, hyp)
# 1. cluster genes into groups based on evenly spaced breaks
# 2. compare the frequency of ``true'' genes to the mid of the cluster

library(plyr)

C <- list(high_mean = matrix(c(0,1,1,0,0,
                               0,-1,1,0,0), 2,5,byrow=T),
          hp_h12    = matrix(c(0,1,1,1,0,
                               0,-1,1,1,0), 2,5,byrow=T),
          lp_h12    = matrix(c(0,1,-1,-1,0,
                               0,-1,-1,-1,0), 2, 5, byrow=T),
          hp_h21    = matrix(c(0,1,1,-1,0,
                               0,-1,1,-1,0), 2,5, byrow=T),
          lp_h21    = matrix(c(0,1,-1,1,0,
                               0,-1,-1,1,0), 2, 5, byrow=T),
          de_p1     = matrix(c(0,1,0,0,0), 1, 5, byrow=T),
          high_or_mid = matrix(c(0,0,1,0,0), 1, 5, byrow=T)
)

binned_probs <- ldply(1:10, function(i){
  
  mcmcfile <- paste("saved_mcmc/mcmc_voom_sim_",i,sep="")
  simfile <- paste("sims/sim_",i,sep="")
  mcmc <- readRDS(mcmcfile)
  sim <- readRDS(simfile)
  
  probs <- mcmc$summaries$probs
  true <- sim$truth$beta
  products <- lapply(C, function(hyp) hyp%*%t(true)>0)
  eval_hyp <- lapply(products, function(product) apply(product, 2, min))
  # str(beta_samp)
  # str(beta_true)
  
  pts <- ldply(1:6, function(j){
    ct <- cut(probs[j,], breaks = 0:10*.1, include.lowest = TRUE)
    out <- sapply(1:10, function(br) mean(eval_hyp[[j]][as.integer(ct)==br]))
    data.frame(hyp=names(C)[j],mids = 0:9*.1 + .05, freq_true = out)
  })
  
  data.frame(pts, sim=i)
  
})

ggplot(binned_probs, aes(x=mids, y=freq_true, group=sim)) +
  geom_line(alpha=.5) + geom_abline(slope=1) +
  facet_wrap(~hyp) +
  theme_bw()
