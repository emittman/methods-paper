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

source("hist-diff.R")
#combine samples of P
s <- readRDS("arrayOutput/chain1_sb_1024.rds")
Ps1 <- s$samples$P
s <- readRDS("arrayOutput/chain2_sb_1024.rds")
Ps2 <- s$samples$P
s <- readRDS("arrayOutput/chain3_sb_1024.rds")
Ps3 <- s$samples$P
s <- readRDS("arrayOutput/chain4_sb_1024.rds")
Ps4 <- s$samples$P
rm(list="s")
library(abind)
P <- abind(Ps1, Ps2, along = 3)
P <- abind(P, Ps3, along=3)
P <- abind(P, Ps4, along=3)
attributes(P)$dimnames[[3]] <- as.character(1:400)
p <- PointwiseMedianHex(P=P, dim1 = 3, dim2 = 4, bins = 30, zmax = 1, alpha=.5)

p1 <- ggplot(df.bnp_est, aes(x=`parental HD`, y=`hybrid`)) + 
  geom_hex(aes(fill=..density..), bins=30)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw()
p2 <- ggplot(df.ind_est, aes(x=`parental HD`, y=`hybrid`)) + 
  geom_hex(aes(fill=..density..), bins=30)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw()
library(cowplot)
no_lgd <- theme_bw(base_size=14) + theme(legend.position = "none")
plot_grid(p2+no_lgd, p[[1]]+no_lgd, p1+no_lgd, nrow=1)

highlow_either.id <- union(union(union(hp12.id, hp21.id),lp12.id),lp21.id)

annot <- annotate(geom = "text", x=0, y=3.8, label=paste(expression(paste(beta[3]," > |",beta[2],"|"))), parse=TRUE)
annot2 <- annotate(geom = "text", x=0, y=-3, label=paste(expression(paste(beta[3]," < -|",beta[2],"|"))), parse=TRUE)

p1 <- ggplot(df.ind_est[highlow_either.id,], aes(x=`parental HD`, y=`hybrid`)) + geom_point(size=.01,alpha=.5) +
  geom_hex(data=df.ind_est, inherit.aes = FALSE, aes(x=`parental HD`, y=`hybrid`, fill=..density..),alpha=.5, bins=32)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw() + geom_abline(slope=c(-1,1), lty=2)+annot + annot2

ylims <- ggplot_build(p1)$layout$panel_ranges[[1]]$y.range
xlims <- ggplot_build(p1)$layout$panel_ranges[[1]]$x.range
p2 <- p[[1]] + geom_point(data=df.bnp_est[highlow_either.id,], inherit.aes = F, 
                    aes(x=`parental HD`, y=`hybrid`), size=.01, alpha=.5) +
  xlim(xlims) +ylim(ylims) + theme_bw(base_size=14) + no_lgd + geom_abline(slope=c(-1,1), lty=2)+ 
  xlab("parental HD") + ylab("hybrid")+annot + annot2
plot_grid(p1+no_lgd, p2, nrow=1) 

p_de <- PointwiseMedianHex(P=s$samples$P, dim1 = 2, dim2 = 3, bins = 30, zmax = 1, alpha=.5)
p3 <- ggplot(df.ind_est[highlow_either.id,], aes(x=`parental HD`, y=`hybrid`)) + geom_point(size=.01) +
  geom_hex(data=df.ind_est, inherit.aes = FALSE, aes(x=`parental HD`, y=`hybrid`, fill=..density..),alpha=.5)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw() + geom_abline(slope=c(-1,1), lty=2)