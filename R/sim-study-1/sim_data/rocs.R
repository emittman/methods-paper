load("../../data/heterosis_design.RData")

library(plyr)
library(dplyr)
library(limma)

grp1 <- 1:3
grp2 <- 4:6
grp3 <- 7:10

roc.df <- function(grp){
  ldply(grp, function(i){
    sim <- readRDS(paste("sim_",i,sep=""))
    mcmc <- readRDS(paste("../saved_mcmc/mcmc_sim_",i,sep=""))
    
    df <- ldply(2:5, function(p){
      ldply(.2*0:5, function(th){
        truth <- abs(sim$truth$beta[,p])>th
        
        #ROC for limma fit
        limma_fit <- lmFit(sim$y, design = X)
        limma_fit <- eBayes(limma_fit)
        tt <- topTreat(fit = limma_fit, coef = p, lfc = th, number=10000)
        tt$g <- row.names(tt)
        # tt <- filter(tt, logFC>0)
        id.limma <- as.numeric(tt$g)
        roc.limma <- data.frame(type="limma",
                                TPR = cumsum(truth[id.limma])/sum(truth),
                                FPR = cumsum(1-truth[id.limma])/sum(1-truth))
        
        #ROC for BNP fit
        p.bnp <- apply(mcmc$samples$beta[p,,],2,function(g){mean(abs(g)>th)})
        id.bnp <- order(p.bnp, decreasing=TRUE)
        roc.bnp <- data.frame(type="bnp",
                              TPR=cumsum((truth)[id.bnp])/sum(truth),
                              FPR = cumsum(1-truth[id.bnp])/sum(1-truth))
        
        data.frame(sim=i, threshold=th, p=p, rbind(roc.bnp, roc.limma)) %>% filter(FPR<.1)
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
  group_by(sim, threshold, p, type, FPR) %>%
  summarize(TPR = max(TPR)) %>%
  mutate(id = paste(type,"_",sim,sep=""))

consolidated$p <- factor(consolidated$p, levels=2:5,
                                  labels=c("parental HD", "hybrid",
                                           "hybrid HD", "flow cell"))
consolidated$type = factor(consolidated$type, levels=c("bnp","limma"),
                                     labels=c("BNP", "limma"))
ssss <- sample(1:10,1)
p1 <- filter(consolidated,threshold<.8, threshold>0, sim==ssss) %>%#,p %in% c("parental HD","hybrid","hybrid HD")) %>%
  # ddply(.(threshold,p,type,FPR), summarise, TPR=mean(TPR)) %>%
  ggplot(aes(FPR,TPR, color=type, linetype=type, group=id)) + geom_line() + 
    facet_grid(threshold~p, scales = "free")+theme_bw(base_size=14)+
  theme(legend.position = c(.16,.78),
        legend.margin = margin(-10,-10,-10,-10,unit="pt"))
p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust=.5))
#AUC
aucs <- arrange(consolidated, sim, threshold, p, type, FPR) %>%
  filter(FPR < .1) %>%
  ddply(.(sim, threshold, p, type), function(x){
    width = diff(x$FPR)
    height = (x$TPR[-1] + x$TPR[-length(x$TPR)])/2
    data.frame(AUC = sum(width*height))
  })
p2 <- filter(aucs, threshold>0,threshold<.8)%>%#, p %in% c("parental HD","hybrid","hybrid HD")) %>%
  ggplot(aes(x=type, y=AUC)) + geom_boxplot(aes(color=type)) +
  geom_line(aes(group=sim), alpha=.5) + facet_grid(threshold~p) +
  theme_bw(base_size=14)+xlab("")+theme(legend.position = "none")
p2
library(cowplot)
plot_grid(p1,p2,ncol=1)
ggsave("../../../figures_tables/ss1-roc-auc.pdf", width=6, height=10)
