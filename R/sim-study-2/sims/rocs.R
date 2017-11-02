load("../../data/heterosis_design.RData")

library(plyr)
library(dplyr)
library(DESeq2)
library(limma)
library(edgeR)

grp1 <- 1:3
grp2 <- 4:6
grp3 <- 7:10

load("../../data/heterosis_design.RData")

roc.df <- function(grp){
  ldply(grp, function(i){
    sim <- readRDS(paste("sim_",i,sep=""))
    mcmc <- readRDS(paste("../saved_mcmc/mcmc_voom_sim_",i,sep=""))
    
    #Fit DESeq2
    deseq.dat <- DESeqDataSetFromMatrix(sim$y, colData = as.data.frame(X[,2:5]),
                                        design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)
    fit.deseq <- DESeq(deseq.dat)
    
    #Fit Voom
    voom.dat <- DGEList(sim$y) %>% calcNormFactors(method="TMM")
    voom.dat.voom <- voom(voom.dat, design = X)
    voom.fit <- lmFit(voom.dat.voom)
    
    # #Fit edgeR
    edgeR.dat <- estimateCommonDisp(voom.dat) %>%
                    estimateTagwiseDisp()
    fit.edgeR <- glmFit(edgeR.dat, design=X)
    
    df <- ldply(2:5, function(p){
      ldply(c(.25,.5,.75), function(th){
        truth <- sim$truth$beta[,p]>th
        
        contr <- rep(0,5)
        contr[p] <- 1
        #ROC for DESeq2
        deseq.tests <- results(fit.deseq, contr, altHypothesis = "greater", lfcThreshold = th)
        id.deseq <- order(deseq.tests$pvalue)
        roc.deseq <- data.frame(type="DESeq2",
                                TPR = cumsum(truth[id.deseq])/sum(truth),
                                FPR = cumsum(1-truth[id.deseq])/sum(1-truth))
        
        #ROC for voom
        voom.fit <- treat(fit = voom.fit, lfc = th)
        voom.top <- topTreat(fit = voom.fit, coef = p, number=10000)
        voom.top$g <- row.names(voom.top)
        voom.top <- filter(voom.top, logFC>0)
        id.voom <- as.numeric(voom.top$g)
        roc.voom <- data.frame(type="voom-limma",
                                TPR = cumsum(truth[id.voom])/sum(truth),
                                FPR = cumsum(1-truth[id.voom])/sum(1-truth))
        
        #ROC for edgeR
        edgeR.tests <- glmTreat(glmfit = fit.edgeR, coef = p, lfc = th)
        edgeR.topG <- topTags(edgeR.tests, n = 10000)[[1]]
        edgeR.topG$g <- row.names(edgeR.topG)
        edgeR.topG <- filter(edgeR.topG, logFC>0)
        id.edgeR <- as.numeric(edgeR.topG$g)
        roc.edgeR <- data.frame(type="edgeR",
                                TPR = cumsum(truth[id.edgeR])/sum(truth),
                                FPR = cumsum(1-truth[id.edgeR])/sum(1-truth))
        
        #ROC for BNP fit
        p.bnp <- apply(mcmc$samples$beta[p,,],2,function(g){mean(g>th)})
        id.bnp <- order(p.bnp, decreasing=TRUE)
        roc.bnp <- data.frame(type="BNP",
                              TPR=cumsum((truth)[id.bnp])/sum(truth),
                              FPR = cumsum(1-truth[id.bnp])/sum(1-truth))
        
        data.frame(sim=i, threshold=th, p=p, rbind(roc.bnp, roc.deseq, roc.voom, roc.edgeR)) %>% filter(FPR<.1)
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
  dplyr::mutate(id = paste(type,"_",sim,sep="")) #%>%
#  dplyr::filter(FPR<.05)

  consolidated %>% #filter(type != "voom-limma") %>%
  # ddply(.(threshold,p,type,FPR), summarise, TPR=mean(TPR)) %>%
  ggplot(aes(FPR,TPR, color=type, linetype=type, group=id)) + 
    geom_line(alpha=.5) + 
    facet_grid(threshold~p, scales = "free")+theme_bw() +
    scale_y_continuous(trans="logit", breaks=c(.05,.25, .5, .75, .9, .95, .99))

ggsave("../../../figures_tables/roc-ss2.pdf",width=7,height=10)
