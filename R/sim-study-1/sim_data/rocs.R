load("../../data/heterosis_design.RData")

library(plyr)
library(dplyr)
library(limma)

roc.df <- ldply(1:10, function(i){
  sim <- readRDS(paste("sim_",i,sep=""))
  mcmc <- readRDS(paste("../saved_mcmc/mcmc_sim_",i,sep=""))
  truth <- sim$truth$beta[,p]>th
  
  df <- ldply(.2*0:5, function(th){
    ldply(2:5, function(p){
      
      #ROC for limma fit
      limma_fit <- lmFit(sim$y, design = X)
      limma_fit <- eBayes(limma_fit)
      tt <- topTreat(fit = limma_fit, coef = p, lfc = th, number=10000)
      tt$g <- row.names(tt)
      tt <- filter(tt, logFC>0)
      id.limma <- as.numeric(tt$g)
      roc.limma <- data.frame(type="limma",
                              TPR = cumsum(truth[id.limma])/sum(truth),
                              FPR = cumsum(1-truth[id.limma])/sum(1-truth))
      
      #ROC for BNP fit
      p.bnp <- apply(mcmc$samples$beta[p,,],2,function(g){mean(g>th)})
      id.bnp <- order(p.bnp, decreasing=TRUE)
      roc.bnp <- data.frame(type="bnp",
                            TPR=cumsum((truth)[id.bnp])/sum(truth),
                            FPR = cumsum(1-truth[id.bnp])/sum(1-truth))
      
      data.frame(threshold=th, p=p, rbind(roc.bnp, roc.limma)) %>% filter(FPR<.1)
    })
  })
  df
})

saveRDS(roc.df, "../data/roc-df.rds")