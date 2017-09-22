source("functions.R")
#demonstrate using generate_log_counts

library(dplyr)
library(limma)
load("data/heterosis_long.RData")
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

sims.ss.1 <- plyr::llply(1:30, function(i){
  y <- generate_log_counts(1000, X, s$samples$P)
  return(y)
})

mse_limma <- plyr::ldply(sims.ss.1, function(sim){
  fit <- limma::lmFit(sim$y, X)
  se_limma <- as.data.frame((coef(fit) - sim$truth$beta)^2)
  colMeans(se_limma)
})

library(ggplot2)
#plot MSE by par over sims
as.data.frame(mse_limma) %>%
  tidyr::gather(key=key, value=value) %>%
  ggplot(aes(x=key,y=value)) + geom_boxplot(aes(fill=key))

#plot ROC curve taking 1-p as ranking of genes
library(pROC)
sim <- sims.ss.1[[1]]
het_pos <- sim$truth$beta[,1]>0
fit <- limma::lmFit(sim$y, X)
fit_eb <- eBayes(fit)
my_roc <- roc(response=het_pos, predictor=fit_eb$p.value[,1], smooth = TRUE)
plot(my_roc)

fits <- lapply(1:30, function(i){
  sim <- sims.ss.1[[i]]
  fit <- limma::lmFit(sim$y, X)
  fit_eb <- eBayes(fit)
  1-fit_eb$p.value[,1]
})

truths <- lapply(1:30, function(i){
  sim <- sims.ss.1[[i]]
  het_pos <- sim$truth$beta[,1]>0
  het_pos
})


plot_rocs <- function(fits, truths, xlimit){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  rocs <- lapply(1:length(fits), function(i){
    fit <- fits[[i]]
    truth <- truths[[i]]
    ord <- order(fit, decreasing=TRUE)
    n_neg <- sum(1-truth)
    n_pos <- sum(truth)
    df <- data.frame(sim=i, true_discovery = cumsum(truth[ord])/n_pos, false_discovery = cumsum(1-truth[ord])/n_neg)
    df
  })
  do.call(rbind, rocs) %>%
    filter(false_discovery<xlimit) %>%
    ggplot(aes(x=false_discovery, y=true_discovery)) + geom_line(aes(group=sim, color=factor(sim)))+
    theme_bw()+geom_abline(slope=1)+
    xlab("FPR")+ylab("TPR")
}

plot_rocs(fits, truths, 1)
