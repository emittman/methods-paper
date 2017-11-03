library(tidyr)
library(plyr)
library(dplyr)

setwd("R/sim-study-1/sim_data/")
load("../../data/heterosis_design.RData")

grp1 <- 1:3
grp2 <- 4:6
grp3 <- 7:10


error.df <- function(grp){
  ldply(grp, function(i){
    sim <- readRDS(paste("sim_",i,sep=""))
    mcmc <- readRDS(paste("../saved_mcmc/mcmc_sim_",i,sep=""))

    #Get BNP fit
    est.df.bnp <- data.frame(t(mcmc$summaries$means_betas))
    names(est.df.bnp) <- c("intercept","parent_hd","hybrid","hybrid_hd","flow_cell")
    est.df.bnp$g <- 1:nrow(est.df.bnp)
    bnp.long <- gather(est.df.bnp, key=par, value=BNP, -g)
    
    #Fit limma
    limma.fit <- lmFit(sim$y, design = X)
    est.df.limma <- data.frame(coef(limma.fit))
    names(est.df.limma)[1] <- "intercept"
    est.df.limma$g <- 1:nrow(est.df.bnp)
    limma.long <- gather(est.df.limma, key=par, value=limma, -g)
    
    #get truth
    truth.df <- data.frame(sim$truth$beta, row.names = NULL)
    names(truth.df)<- c("intercept","parent_hd","hybrid","hybrid_hd","flow_cell")
    truth.df$g <- 1:nrow(truth.df)
    truth.long <- gather(truth.df, key=par, value=truth, -g)
    
    merged.all <- merge(bnp.long, truth.long, by=c("g","par"))
    merged.all <- merge(merged.all, limma.long, by=c("g","par"))
    
    means.df <- gather(merged.all, key=type, value=est, BNP, limma)
    
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
  dplyr::mutate(par = factor(par, levels=c("intercept","parent_hd","hybrid","hybrid_hd","flow_cell"),
                             labels=c("intercept","parental HD", "hybrid", "hybrid HD", "flow cell"))) %>%
  gather(key=error_type, value=value, MSPE, MAPE)
  
combined$type <- factor(combined$type, levels=c("BNP","limma"))

combined %>%
  dplyr::filter(par != "intercept", error_type=="MAPE") %>%
  ggplot(aes(type, y=value, color=type)) + 
  geom_boxplot(aes(fill=type), alpha=.5) + 
  facet_wrap(~par, scales="free")+
  theme_bw(base_size=14) + ylab("") + xlab("") + ggtitle("Mean Absolute Prediction Error")

p3 <- combined %>%
  dplyr::filter(par != "intercept", error_type=="MSPE", type != "voom-limma") %>%
  ggplot(aes(type, y=value)) +
  geom_line(aes(group=sim))+
  geom_boxplot(aes(color=type), alpha=.5) + 
  facet_wrap(~par, scales="free")+
  theme_bw(base_size=14) + ylab("") + xlab("") + ggtitle("Mean Squared Prediction Error")+
  theme(legend.position = "none")

ggsave("../../../figures_tables/ss1-mspe.pdf", width=6, height=5)
