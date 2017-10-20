library(limma)

load("R/data/heterosis_design.RData")

library(edgeR)
library(dplyr)
library(DESeq2)
deseq.dat <- DESeqDataSetFromMatrix(sim$y, colData = as.data.frame(X[,2:5]),
                                    design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)

fit.deseq <- DESeq(deseq.dat)
df.deseq2 <- as.data.frame(coef(fit.deseq))
names(df.deseq2)[1] <- "intercept"

#This shows the priors being used on betas by DESeq2
mles <- estimateMLEForBetaPriorVar(fit.deseq)
priors <- estimateBetaPriorVar(mles, betaPriorMethod = "quantile")

dge <- DGEList(sim$y) %>% calcNormFactors(method="TMM")

voom_dge <- voom(dge, design = X)
fit <- lmFit(voom_dge)
df.voom <- data.frame(coef(fit))
names(df.voom)[1] <- "intercept"
fit.eb <- eBayes(fit)

dge <- estimateCommonDisp(dge) %>%
  estimateTagwiseDisp()
fit.edgeR <- glmFit(dge, design=X)

df.edgeR <- data.frame(coef(fit.edgeR))
names(df.edgeR)[1] <- "intercept"


#rocs for beta3
deseq.tests <- results(fit.deseq, c(0,0,1,0,0), altHypothesis = "greater")
deseq.topG <- order(deseq.tests$pvalue)

edgeR.tests <- glmLRT(glmfit = fit.edgeR, coef = 3)
edgeR.topG <- topTags(edgeR.tests, n = 10000)[[1]]
edgeR.topG <- edgeR.topG[which(edgeR.topG$logFC>0),]
edgeR.topG <- as.numeric(row.names(edgeR.topG))


voom.top <- topTable(fit.eb, coef = 3, number = 10000, confint = TRUE)
voom.top <- voom.top[which(voom.top$logFC>0),]
voom.topG <- as.numeric(row.names(voom.top))

bnp.topG <- order(s$summaries$probs[7,], decreasing = TRUE)
# which.pos <- which(s$summaries$means_betas[3,bnp.topG]>-0.1)
# bnp.topG <- bnp.topG[which.pos]

# bnp2.topG <- order(s2$summaries$probs[7,], decreasing = TRUE)
# which.pos2 <- which(s2$summaries$means_betas[3,bnp2.topG]>-0.1)
# bnp2.topG <- bnp2.topG[which.pos2]
gene_id <- which(ses.hyb<.007)
bnp2.topG <- bnp.topG[which(bnp.topG %in% gene_id)]


true.topG <- sim$truth$beta[,3]>0
total.true <- sum(true.topG)

edgeR.tpr <- cumsum(true.topG[edgeR.topG])/total.true
edgeR.fpr <- cumsum(1-true.topG[edgeR.topG])/(10000-total.true)
edgeR.roc <- data.frame(fpr = edgeR.fpr,
                        tpr = edgeR.tpr,
                        method="edgeR")

deseq2.tpr <- cumsum(true.topG[deseq.topG])/total.true
deseq2.fpr <- cumsum(1-true.topG[deseq.topG])/(10000-total.true)
deseq2.roc <- data.frame(fpr = deseq2.fpr,
                        tpr = deseq2.tpr,
                        method="deseq2")


voom.tpr <- cumsum(true.topG[voom.topG])/total.true
voom.fpr <- cumsum(1-true.topG[voom.topG])/(10000-total.true)
voom.roc <- data.frame(fpr = voom.fpr,
                       tpr = voom.tpr,
                       method="voom")

bnp.tpr <- cumsum(true.topG[bnp.topG])/total.true
bnp.fpr <- cumsum(1-true.topG[bnp.topG])/(10000-total.true)
bnp.roc <- data.frame(fpr = bnp.fpr,
                      tpr = bnp.tpr,
                      method="bnp")

bnp2.tpr <- cumsum(true.topG[bnp2.topG])/total.true
bnp2.fpr <- cumsum(1-true.topG[bnp2.topG])/(10000-total.true)
bnp2.roc <- data.frame(fpr = bnp2.fpr,
                      tpr = bnp2.tpr,
                      method="bnp_low_se")

p1 <- rbind(bnp.roc,voom.roc,deseq2.roc,edgeR.roc,bnp2.roc) %>%
  filter(fpr<.05) %>%
  ggplot(aes(fpr,tpr, color=method,group=method)) + geom_step(size=1) +
  geom_abline(slope=1) + ggtitle("ROCs for 'hybrid effect'")

#rocs for de
deseq.tests <- results(fit.deseq, c(0,1,0,0,0), altHypothesis = "greater")
deseq.topG <- order(deseq.tests$pvalue)

edgeR.tests <- glmLRT(glmfit = fit.edgeR, coef = 2)
edgeR.topG <- topTags(edgeR.tests, n = 10000)[[1]]
edgeR.topG <- edgeR.topG[which(edgeR.topG$logFC>0),]
edgeR.topG <- as.numeric(row.names(edgeR.topG))

voom.top <- topTable(fit.eb, coef = 2, number = 10000, confint = TRUE)
voom.top <- voom.top[which(voom.top$logFC>0),]
voom.topG <- as.numeric(row.names(voom.top))

bnp.topG <- order(s$summaries$probs[6,], decreasing = TRUE)
gene_id <- which(ses < .008)
bnp2.topG <- bnp.topG[bnp.topG %in% gene_id]
# bnp.topG <- bnp.topG[which(s$summaries$means_betas[2,bnp.topG]>-.1)]

# bnp2.topG <- order(s2$summaries$probs[6,], decreasing = TRUE)
# bnp2.topG <- bnp2.topG[which(s2$summaries$means_betas[2,bnp2.topG]>-.1)]

true.topG <- sim$truth$beta[,2]>0
total.true <- sum(true.topG)

edgeR.tpr <- cumsum(true.topG[edgeR.topG])/total.true
edgeR.fpr <- cumsum(1-true.topG[edgeR.topG])/(10000-total.true)
edgeR.roc <- data.frame(fpr = edgeR.fpr,
                        tpr = edgeR.tpr,
                        method="edgeR")

deseq2.tpr <- cumsum(true.topG[deseq.topG])/total.true
deseq2.fpr <- cumsum(1-true.topG[deseq.topG])/(10000-total.true)
deseq2.roc <- data.frame(fpr = deseq2.fpr,
                         tpr = deseq2.tpr,
                         method="deseq2")

voom.tpr <- cumsum(true.topG[voom.topG])/total.true
voom.fpr <- cumsum(1-true.topG[voom.topG])/(10000-total.true)
voom.roc <- data.frame(fpr = voom.fpr,
                       tpr = voom.tpr,
                       method="voom")

bnp.tpr <- cumsum(true.topG[bnp.topG])/total.true
bnp.fpr <- cumsum(1-true.topG[bnp.topG])/(10000-total.true)
bnp.roc <- data.frame(fpr = bnp.fpr,
                      tpr = bnp.tpr,
                      method="bnp")

bnp2.tpr <- cumsum(true.topG[bnp2.topG])/total.true
bnp2.fpr <- cumsum(1-true.topG[bnp2.topG])/(10000-total.true)
bnp2.roc <- data.frame(fpr = bnp2.fpr,
                      tpr = bnp2.tpr,
                      method="bnp (low std.errs)")

p2 <- rbind(bnp.roc,voom.roc,deseq2.roc,edgeR.roc,bnp2.roc) %>%
  filter(fpr<.05) %>%
  ggplot(aes(fpr,tpr, color=method,group=method)) + geom_step(size=1) +
  geom_abline(slope=1)+
  ggtitle("ROCs for B73 > Mo17")

library(cowplot)
lg <- get_legend(p1)
no.lg <- theme(legend.position = "none")
plot_grid(p1+no.lg,p2+no.lg,lg, nrow=1, rel_widths = c(1,1,.2))


df.bnp <- data.frame(t(s$summaries$means_betas))
names(df.bnp) <- names(df.voom)

df.bnp2 <- data.frame(t(s2$summaries$means_betas))
names(df.bnp2) <- names(df.voom)

df.voom$method <- "voom"
df.bnp$method <- "bnp"
df.bnp2$method <- "bnp_novoom"
df.edgeR$method <- "edgeR"
df.deseq2$method <- "DESeq2"

names(df.edgeR)[1:5] <- names(df.voom)[1:5]

df.truth <- data.frame(sim$truth$beta)
df.truth$method <- 'truth'
names(df.truth)[1:5] <- names(df.voom)[1:5]
library(ggplot2)
library(tidyr)
rbind(df.voom,df.bnp,df.edgeR, df.truth) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=parent_hd, y=hybrid, color=method)) + geom_point(alpha=.01) +
  facet_wrap(~method) + 
  theme_bw()+ stat_function(fun=abs)

rbind(df.voom,df.bnp,df.bnp2) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=intercept, y=hybrid, color=method)) + geom_density_2d(alpha=.2, n=100, binwidth=.01) +
  facet_wrap(~method, scales="free") + theme_bw()

rbind(df.voom,df.bnp,df.truth) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=parent_hd, y=flow_cell, color=method)) + geom_density_2d(alpha=.5, n=50, binwidth=.04) +
  facet_wrap(~method) + theme_bw()

rbind(df.voom,df.bnp) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=intercept, y=parent_hd, color=method)) + geom_density_2d(alpha=.5, n=50, binwidth=.01) +
  facet_wrap(~method) + theme_bw()

rbind(df.voom,df.bnp,df.edgeR,df.truth) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=intercept, y=hybrid_hd, color=method)) + geom_density_2d(alpha=.5, n=50, binwidth=.04) +
  facet_wrap(~method) + theme_bw()

rbind(df.voom,df.bnp,df.edgeR,df.truth) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=flow_cell, y=hybrid_hd, color=method)) + geom_density_2d(alpha=.5, n=100, binwidth=.5) +
  facet_wrap(~method) + theme_bw()

id.tails <- 1:10000#which(abs(sim$truth$beta[,2])>1.25)

mse.voom <- sapply(2:5, function(p){
  sq.err <- (df.voom[id.tails,p] - df.truth[id.tails,p])^2
  mean(sq.err)
})
mse.bnp <- sapply(2:5, function(p){
  sq.err <- (df.bnp[id.tails,p] - df.truth[id.tails,p])^2
  mean(sq.err)
})
mse.bnp2 <- sapply(2:5, function(p){
  sq.err <- (df.bnp2[id.tails,p] - df.truth[id.tails,p])^2
  mean(sq.err)
})
mse.edger <- sapply(2:5, function(p){
  sq.err <- (df.edgeR[id.tails,p] - df.truth[id.tails,p])^2
  mean(sq.err)
})
mse.deseq2 <- sapply(2:5, function(p){
  sq.err <- (df.deseq2[id.tails,p] - df.truth[id.tails,p])^2
  mean(sq.err)
})

data.frame(method=c("voom","bnp","bnp2","edgeR","DESeq2"),
           rbind(mse.voom,mse.bnp,mse.bnp2,mse.edger,mse.deseq2))

df2.voom <- rbind(df.voom, df.truth) %>%
  mutate(gene = rep(1:10000,2)) %>%
  gather(key=par, value=value, -method, -gene) %>%
  spread(key=method, value=value)

df2.bnp <- rbind(df.bnp, df.truth) %>%
  mutate(gene = rep(1:10000,2)) %>%
  gather(key=par, value=value, -method, -gene) %>%
  spread(key=method, value=value)

df2.deseq2 <- rbind(df.deseq2, df.truth) %>%
  mutate(gene = rep(1:10000,2)) %>%
  gather(key=par, value=value, -method, -gene) %>%
  spread(key=method, value=value)

df2.edgeR <- rbind(df.edgeR, df.truth) %>%
  mutate(gene = rep(1:10000,2)) %>%
  gather(key=par, value=value, -method, -gene) %>%
  spread(key=method, value=value)

df2 <- merge(df2.deseq2, df2.edgeR, 
             by=c("gene","par","truth"))
df2 <- merge(df2, df2.bnp, by=c("gene","par","truth"))
df2 <- merge(df2, df2.voom, by=c("gene","par","truth"))

df2 <- gather(df2, key=method, value=value, 4:7)

filter(df2, par != "intercept") %>%
ggplot(aes(truth, value, color=method, group=method)) +
  geom_point(alpha=.2, pch=".") + facet_wrap(~par, scales="free") + 
  geom_abline(slope=1) +
  geom_smooth()

spread(df2, key=method, value=value) %>%
  ggplot(aes(voom, edgeR)) + geom_hex() +
  scale_fill_continuous(trans="log", low="white", high="black") +
  theme_bw() + facet_wrap(~par, scales="free") + geom_abline(slope=1)

#plot draw of 'P' vs truth (as contour plot)
ggplot(data=NULL, aes(x = s$samples$P[,3,27], y = s$samples$P[,4,27], weights = s$samples$P[,1,27]))+
  geom_hex(bins=20) + scale_fill_continuous(trans="log", low="white", high="black") +
  theme_bw() + 
  geom_density_2d(data=df.truth, inherit.aes=FALSE, aes(x=parent_hd, y=hybrid), n=100, binwidth=.03)

df.P <- with(s$samples, rbind(data.frame(parent_hd = P[,3,30],
                                         hybrid=P[,4,30],
                                         weight=P[,1,30],
                                         facet='1'),
                              data.frame(parent_hd = P[,3,40],
                                         hybrid=P[,4,40],
                                         weight=P[,1,40],
                                         facet='2'),
                              data.frame(parent_hd = P[,3,50],
                                         hybrid=P[,4,50],
                                         weight=P[,1,50],
                                         facet='3'),
                              cbind(df.truth[,c('parent_hd','hybrid')],
                                    data.frame(weight=rep(.0001,10000),
                                               facet='truth'))
             ))

ggplot(df.P, aes(parent_hd,hybrid,weight=weight)) +
  geom_hex(bins=15) + scale_fill_continuous(trans="log",
                                     low="white", high="darkblue",
                                     breaks=c(10^(-7:-1))) +
  theme_bw() + facet_wrap(~facet)

# plots of 'P'
P <- s$samples$P
pl1 <- PointwiseMedianHex(P, 2, 3, 20)
pl2 <- PointwiseMedianHex(P, 2, 4, 20)
pl3 <- PointwiseMedianHex(P, 2, 5, 20)
pl4 <- PointwiseMedianHex(P, 2, 6, 20)
pl5 <- PointwiseMedianHex(P, 3, 4, 20)
pl6 <- PointwiseMedianHex(P, 3, 5, 20)
pl7 <- PointwiseMedianHex(P, 3, 6, 20)
pl8 <- PointwiseMedianHex(P, 4, 5, 20)
pl9 <- PointwiseMedianHex(P, 4, 6, 20)
P[,7,] <- sqrt(1/P[,7,])
pl10 <- PointwiseMedianHex(P, 2, 7, 20)

library(cowplot)
plot_grid(pl1[[1]]+nolgd,pl1[[2]]+nolgd,
          pl2[[1]]+nolgd,pl2[[2]]+nolgd,
          pl3[[1]]+nolgd,pl3[[2]]+nolgd, ncol=2)
plot_grid(pl4[[1]]+nolgd,pl4[[2]]+nolgd,
          pl5[[1]]+nolgd,pl5[[2]]+nolgd,
          pl6[[1]]+nolgd,pl6[[2]]+nolgd, ncol=2)
plot_grid(pl7[[1]]+nolgd,pl7[[2]]+nolgd,
          pl8[[1]]+nolgd,pl8[[2]]+nolgd,
          pl9[[1]]+nolgd,pl9[[2]]+nolgd,
          pl10[[1]]+nolgd,pl10[[2]]+nolgd,ncol=2)

# binned probabilities
sim <- readRDS("sims/sim_1")
binned_probs <- data.frame(bins = cut(s$summaries$probs[6,],
                                      breaks = 0:20*.05,
                                      include.lowest = TRUE))
binned_probs$mid <- c(0:19*.05 + .025)[match(binned_probs$bins, sort(unique(binned_probs$bins)))]
binned_probs$true <- sim$truth$beta[,3]>0
ddply(binned_probs, .(mid), summarize, true_proportion=mean(true), n=length(true)) %>%
  ggplot(aes(mid, true_proportion,fill=n)) + geom_bar(stat = "identity")+
  geom_abline(slope=1)

binned_probs <- data.frame(bins = cut(s$summaries$probs[7,],
                               breaks = 0:20*.05,
                               include.lowest = TRUE))
binned_probs$mid <- c(0:19*.05 + .025)[match(binned_probs$bins, sort(unique(binned_probs$bins)))]
binned_probs$true <- sim$truth$beta[,2]>0
ddply(binned_probs, .(mid), summarize, true_proportion=mean(true), n=length(true)) %>%
  ggplot(aes(mid, true_proportion,fill=n)) + geom_bar(stat = "identity")+
  geom_abline(slope=1)

binned_probs <- data.frame(bins = cut(s$summaries$probs[3,],
                                      breaks = 0:20*.05,
                                      include.lowest = TRUE))
binned_probs$mid <- c(0:19*.05 + .025)[match(binned_probs$bins, sort(unique(binned_probs$bins)))]
binned_probs$true <- sim$truth$beta[,4]>0
ddply(binned_probs, .(mid), summarize, true_proportion=mean(true), n=length(true)) %>%
  ggplot(aes(mid, true_proportion,fill=n)) + geom_bar(stat = "identity")+
  geom_abline(slope=1)
