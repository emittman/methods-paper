#Comparison of estimates and posterior expectation of \mathcal{P}

load("data/ind-est-heterosis.RData")

str(ind_est)

df.ind_est <- with(ind_est, data.frame(t(beta), sigma=sqrt(sigma2)))
names(df.ind_est) <- c("intercept", "parental HD", "hybrid", "hybrid HD", "flow cell", "sigma")
library(ggplot2)
library(GGally)
library(dplyr)
my_hex <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_hex(aes(fill=..density..), ..., bins=30, color = NA)+
    scale_fill_continuous(trans="log",low="white",high="darkviolet")+
    theme_bw()
}

select(df.ind_est, -sigma) %>%
  ggpairs(lower=list(continuous=my_hex))

load("paschold_analysis/combined_summaries.RData")
df.bnp_est <- with(summaries, data.frame(t(means_betas), sigma=means_sigmas))
names(df.bnp_est) <- c("intercept", "parental HD", "hybrid", "hybrid HD", "flow cell", "sigma")

select(df.bnp_est, -sigma) %>%
  ggpairs(lower=list(continuous=  my_hex))

library(DESeq2)
load("data/heterosis_counts.RData")
load("data/heterosis_design.RData")

colnames(my_dat_wide) <- NULL
deseq.dat <- DESeqDataSetFromMatrix(as.matrix(my_dat_wide[,-1]), colData = as.data.frame(X[,-1]),
                                    design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)
fit.deseq <- DESeq(deseq.dat, betaPrior = TRUE)
df.deseq_est <- as.data.frame(coef(fit.deseq))
names(df.deseq_est) <- names(df.bnp_est)[1:5]
fit.deseq2 <- DESeq(deseq.dat, betaPrior = FALSE)
df.deseq_est2 <- as.data.frame(coef(fit.deseq2))
names(df.deseq_est2) <- names(df.bnp_est)[1:5]
ggpairs(data=df.deseq_est, lower=list(continuous=my_hex))

nolg <- theme(legend.position = "none")

library(cowplot)

df.comb <- rbind(#cbind(df.ind_est[,-6], method="voom weights"),
                 cbind(df.bnp_est[,-6], method="BNP"),
                 cbind(df.deseq_est2, method="DESeq2 (unshrunk)"),
                 cbind(df.deseq_est, method="DESeq2 (shrunk)"))
p1 <- ggplot(data=df.comb, aes(x=`parental HD`, y=hybrid))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p2 <- ggplot(data=df.comb, aes(x=`parental HD`, y=`flow cell`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p3 <- ggplot(data=df.comb, aes(x=`parental HD`, y=`hybrid HD`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p3 <- ggplot(data=df.comb, aes(x=`parental HD`, y=`hybrid HD`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p4 <- ggplot(data=df.comb, aes(x=`hybrid`, y=`hybrid HD`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p5 <- ggplot(data=df.comb, aes(x=`hybrid`, y=`flow cell`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

p6 <- ggplot(data=df.comb, aes(x=`hybrid HD`, y=`flow cell`))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="gray",high="midnightblue")+
  theme_bw(base_size=12)+facet_wrap(~method, scales="fixed") + nolg

plot_grid(p1,p3,p2,p4,p5,p6, nrow=3)

xlims <- ggplot_build(p1)$layout$panel_ranges[[1]]$x.range
ylims <- ggplot_build(p1)$layout$panel_ranges[[1]]$y.range

ggb <- ggplot_build(p1)
str(ggb$data)
p2 <- ggplot(data=df.bnp_est, aes(x=`parental HD`, y=hybrid))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw() + xlim(xlims) + ylim(ylims)
p3 <- ggplot(data=df.deseq_est, aes(x=`parental HD`, y=hybrid))+
  geom_hex(aes(fill=..density..), bins=30, color = NA)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw()+ xlim(xlims) + ylim(ylims)
