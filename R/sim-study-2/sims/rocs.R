load("../../data/heterosis_design.RData")

library(plyr)
library(dplyr)
library(DESeq2)

grp1 <- 1:3
grp2 <- 4:6
grp3 <- 7:10

load("../../data/heterosis_design.RData")

roc.df <- function(grp){
  ldply(grp, function(i){
    sim <- readRDS(paste("sim_",i,sep=""))
    mcmc <- readRDS(paste("../saved_mcmc/mcmc_voom_sim_",i,sep=""))
    deseq.dat <- DESeqDataSetFromMatrix(sim$y, colData = as.data.frame(X[,2:5]),
                                        design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)
    fit.deseq <- DESeq(deseq.dat)

    df <- ldply(2:5, function(p){
      ldply(c(.25,.5,1), function(th){
        truth <- sim$truth$beta[,p]>th
        
        contr <- rep(0,5)
        contr[p] <- 1
        #ROC for DESeq2
        deseq.tests <- results(fit.deseq, contr, altHypothesis = "greater", lfcThreshold = th)
        id.deseq <- order(deseq.tests$pvalue)
        roc.deseq <- data.frame(type="DESeq2",
                                TPR = cumsum(truth[id.deseq])/sum(truth),
                                FPR = cumsum(1-truth[id.deseq])/sum(1-truth))
        
        # deseq2_fit <- lmFit(sim$y, design = X)
        # limma_fit <- eBayes(limma_fit)
        # tt <- topTreat(fit = limma_fit, coef = p, lfc = th, number=10000)
        # tt$g <- row.names(tt)
        # tt <- filter(tt, logFC>0)
        # id.limma <- as.numeric(tt$g)
        
        
        #ROC for BNP fit
        p.bnp <- apply(mcmc$samples$beta[p,,],2,function(g){mean(g>th)})
        id.bnp <- order(p.bnp, decreasing=TRUE)
        roc.bnp <- data.frame(type="bnp",
                              TPR=cumsum((truth)[id.bnp])/sum(truth),
                              FPR = cumsum(1-truth[id.bnp])/sum(1-truth))
        
        data.frame(sim=i, threshold=th, p=p, rbind(roc.bnp, roc.deseq)) %>% filter(FPR<.1)
      })
    })
    df
  })
}
out1 <- roc.df(grp1)
saveRDS(out1, "../data/roc-df1.rds")
out2 <- roc.df(grp2)
saveRDS(out2, "../data/roc-df2.rds")
out3 <- roc.df(grp3)
saveRDS(out3, "../data/roc-df3.rds")

out1 <- readRDS("roc-df1.rds")
out2 <- readRDS("roc-df2.rds")
out3 <- readRDS("roc-df3.rds")

library(ggplot2)
consolidated <- rbind(out1,out2,out3) %>%
  dplyr::group_by(sim, threshold, p, type, FPR) %>%
  dplyr::summarize(TPR = max(TPR)) %>%
  dplyr::mutate(id = paste(type,"_",sim,sep="")) %>%
  dplyr::filter(FPR<.05)

  consolidated %>%
  # ddply(.(threshold,p,type,FPR), summarise, TPR=mean(TPR)) %>%
  ggplot(aes(FPR,TPR, color=type, linetype=type, group=id)) + geom_line(alpha=.5) + 
    facet_grid(threshold~p, scales = "free")+theme_bw()
