#truth of hypotheses
setwd("sim_data/")
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
          de_p1     = matrix(c(0,1,0,0,0), 1, 5, byrow=T)
)

library(plyr)

all_rocs <- ldply(1:10, function(j){
  sim <- readRDS(paste("sim_", j, sep=""))
  true <- sim$truth$beta
  products <- lapply(C, function(hyp) hyp%*%t(true)>0)
  eval_hyp <- lapply(products, function(product) apply(product, 2, min))
  
  mcmc <- readRDS(paste("../saved_mcmc/mcmc_sim_",j,sep=""))
  
  rocs <- ldply(1:6, function(i){
    ord_hm <- order(mcmc$summaries$probs[i,], decreasing =TRUE)
    true_pos_rate <- cumsum(eval_hyp[[i]][ord_hm])/sum(eval_hyp[[i]])
    false_pos_rate <- cumsum(1-eval_hyp[[i]][ord_hm])/sum(1-eval_hyp[[i]])
    data.frame(sim = j,
               hypothesis = names(C)[[i]],
               TPR = true_pos_rate,
               FPR = false_pos_rate)
  })
  rocs
})

save(all_rocs, file="../summaries/all_rocs.RData")
load("../summaries/all_rocs.RData")


library(limma)
y1 <- readRDS("sim_data/sim_1")$y
load("../data/heterosis_long.RData")
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

limma_fit <- lmFit(y1, design = X)
limma_fit <- eBayes(limma_fit)
tt <- topTable(fit = limma_fit, number = 10000, coef = 2)
tt$g <- row.names(tt)
tt <- filter(tt, logFC>0)
id_de_1 <- as.numeric(tt$g)
roc_de_1 <- data.frame(type="limma",
                       TPR = cumsum(eval_hyp[[6]][id_de_1])/sum(eval_hyp[[6]]),
                       FPR = cumsum(1-eval_hyp[[6]][id_de_1])/sum(1-eval_hyp[[6]]))

roc_de_bnp <- filter(all_rocs, sim==1, hypothesis=="de_p1") %>%
  select(TPR,FPR) %>% mutate(type="bnp")

rbind(roc_de_1,roc_de_bnp) %>% filter(FPR<.1)%>%
  ggplot(aes(FPR, TPR, color=type)) + geom_line() + geom_abline(slope=1, lty=2)



library(ggplot2)
library(dplyr)
filter(all_rocs, FPR<.25, sim<8) %>%
ggplot(aes(FPR, TPR)) + geom_line() + geom_abline(slope=1, lty=2) +
  facet_grid(sim~hypothesis)


sim1 <- readRDS("sim_1")
true1 <- sim1$truth$beta
products <- lapply(C, function(hyp) hyp%*%t(true1)>0)
eval_hyp <- lapply(products, function(product) apply(product, 2, min))

mcmc1 <- readRDS("../saved_mcmc/mcmc_sim_1")
ord_hm <- order(mcmc1$summaries$probs[6,], decreasing =FALSE)

true_pos_rate <- cumsum((1-eval_hyp[[6]])[ord_hm])/sum(1-eval_hyp[[6]])
false_pos_rate <- cumsum(eval_hyp[[6]][ord_hm])/sum(eval_hyp[[6]])

plot(false_pos_rate, true_pos_rate)
abline(0,1, col="red")
