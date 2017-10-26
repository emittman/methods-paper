library(tidyr)
library(plyr)
library(dplyr)
library(DESeq2)

setwd("R/sim-study-2/sims/")
load("../../data/heterosis_design.RData")

grp1 <- 1:3
grp2 <- 4:6
grp3 <- 7:10

load("../../data/heterosis_design.RData")

error.df <- function(grp){
  ldply(grp, function(i){
    sim <- readRDS(paste("sim_",i,sep=""))
    mcmc <- readRDS(paste("../saved_mcmc/mcmc_voom_sim_",i,sep=""))

    est.df.bnp <- data.frame(t(mcmc$summaries$means_betas))
    names(est.df.bnp) <- c("intercept","parent_hd","hybrid","hybrid_hd","flow_cell")
    est.df.bnp$g <- 1:nrow(est.df.bnp)
    bnp.long <- gather(est.df.bnp, key=par, value=BNP, -g)
    
    deseq.dat <- DESeqDataSetFromMatrix(sim$y, colData = as.data.frame(X[,2:5]),
                                        design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)
    
    fit.deseq <- DESeq(deseq.dat)
    est.df.deseq2 <- as.data.frame(coef(fit.deseq))
    names(est.df.deseq2)[1] <- "intercept"
    est.df.deseq2$g <- 1:nrow(est.df.bnp)
    ds2.long <- gather(est.df.deseq2, key=par, value=DESeq2, -g)
    
    #Fit Voom
    voom.dat <- DGEList(sim$y) %>% calcNormFactors(method="TMM")
    voom.dat.voom <- voom(voom.dat, design = X)
    voom.fit <- lmFit(voom.dat.voom)
    est.df.voom <- data.frame(coef(voom.fit))
    names(est.df.voom)[1] <- "intercept"
    est.df.voom$g <- 1:nrow(est.df.bnp)
    voom.long <- gather(est.df.voom, key=par, value=voom, -g)
    
    truth.df <- data.frame(sim$truth$beta)
    names(truth.df)[1] <- "intercept"
    truth.df$g <- 1:nrow(truth.df)
    truth.long <- gather(truth.df, key=par, value=truth, -g)
    
    merged.bnp <- merge(bnp.long, truth.long, by=c("g","par"))
    merged.bnp <- merge(merged.bnp, voom.long, by=c("g","par"))
    merged.all <- merge(merged.bnp, ds2.long, by=c("g","par"))
    
    means.df <- gather(merged.all, key=type, value=est, BNP, DESeq2, voom)
    
    out <- means.df %>% 
      dplyr::group_by(g, par, type) %>%
      dplyr::summarize(sq.error = (est-truth)^2,
                abs.error = abs(est-truth)) %>%
      dplyr::group_by(par, type) %>%
      dplyr::summarize(MSPE = mean(sq.error),
                MAPE = mean(abs.error)) %>%
      dplyr::mutate(sim=i)
    
    return(out)
  })
}

out1 <- error.df(grp1)
saveRDS(out1, "../data/error-df1.rds")
out2 <- error.df(grp2)
saveRDS(out2, "../data/error-df2.rds")
out3 <- error.df(grp3)
saveRDS(out3, "../data/error-df3.rds")

out1 <- readRDS("../data/error-df1.rds")
out2 <- readRDS("../data/error-df2.rds")
out3 <- readRDS("../data/error-df3.rds")

# out2 <- means.df %>% 
#   dplyr::group_by(g, par, type) %>%
#   dplyr::summarize(sq.error = (est-truth)^2,
#                    abs.error = abs(est-truth))
# 
# ggplot(filter(out2, par != "intercept"), aes(x=sq.error, group=type)) + geom_density(aes(fill=type), alpha=.2, n=1024) +
#   facet_wrap(~par, scales="free") +
#   scale_x_continuous(trans="log", breaks=c(.00001,.0001,.001,.01,.1))

combined <- rbind(out1,out2,out3) %>%
  dplyr::mutate(par = factor(par, levels=c("intercept","parent_hd","hybrid","hybrid_hd","flow_cell"))) %>%
  gather(key=error_type, value=value, MSPE, MAPE)
  
combined %>%
  dplyr::filter(par != "intercept", error_type=="MAPE") %>%
  ggplot(aes(factor(type), y=value)) + 
  geom_boxplot() + 
  facet_wrap(~par, scales="free")+
  theme_bw(base_size=14) + ylab("") + xlab("") + ggtitle("Mean Absolute Prediction Error")

combined %>%
  dplyr::filter(par != "intercept", error_type=="MSPE") %>%
  ggplot(aes(factor(type), y=value)) + 
  geom_boxplot() + 
  facet_wrap(~par, scales="free")+
  theme_bw(base_size=14) + ylab("") + xlab("") + ggtitle("Mean Squared Prediction Error")

ggsave("../../../figures_tables/ss2-mspe.pdf", width=6, height=5)
