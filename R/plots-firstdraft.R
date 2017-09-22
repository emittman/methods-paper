source("../../chapter2/R/util/rhat.R")

s1 <- readRDS("arrayOutput/chain1_sb.rds")
s2 <- readRDS("arrayOutput/chain2_sb.rds")
s3 <- readRDS("arrayOutput/chain3_sb.rds")
s4 <- readRDS("arrayOutput/chain4_sb.rds")
class(s1) <- "myMcmcObj"
class(s2) <- "myMcmcObj"
class(s3) <- "myMcmcObj"
class(s4) <- "myMcmcObj"

niter <- length(s1$samples$alpha)

rhat <- gelman.factors(s1,s2,s3,s4, n_iter=niter)

rhat.df <- data.frame(t(rhat[[1]]), rhat[[2]])
names(rhat.df) <- c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","sigma")
gather(rhat.df, key=par,value=value) %>%
  ggplot(aes(x=value)) + geom_histogram(bins=100) + facet_wrap(~par,scales="free") +
  theme_classic()

comb <- combine_chains_summaries(s1, s2, s3, s4)

load("data/ind-est-heterosis.RData")

library(ggplot2)

G <- dim(comb$means_betas)[2]

bnp <- data.frame(t(comb$means_betas)) %>%
  mutate(sigma=comb$means_sigmas,
         type="bnp",
         gene.id=1:G)
ols <- data.frame(t(ind_est$beta)) %>%
  mutate(sigma=sqrt(ind_est$sigma2),
         type="ols",
         gene.id=1:G)

names(bnp)[1:5] <- c("intercept","par_hd","hyb","hyb_hd","flowcell")
names(ols)[1:5] <- c("intercept","par_hd","hyb","hyb_hd","flowcell")

library(dplyr)
library(tidyr)
plot.df <- rbind(bnp,ols) %>%
  gather(key=par, value=value, 1:6) %>%
  spread(key=type, value=value)

sam <- head(plot.df)

ggplot(plot.df, aes(bnp,ols)) + geom_hex() + facet_wrap(~par)+
  geom_abline(slope=1) + theme_classic()+
  scale_fill_continuous(low="white",high="black", trans="log",
                        breaks=c(1,10,100,1000,10000))

plot.df2 <- plot.df %>%
  gather(key=type,value=value, 3:4) %>%
  spread(key=par, value=value)

is <- ggplot(plot.df2, aes(intercept,sigma)) + geom_hex() + facet_grid(.~type,scales="fixed")+
  geom_smooth(formula=y~s(x,k=20),method="gam") + theme_classic()+
  scale_fill_continuous(low="white",high="black", trans="log",
                        breaks=c(1,10,100,1000,10000))

hs <- ggplot(plot.df2, aes(hyb,sigma)) + geom_hex() + facet_grid(.~type,scales="fixed")+
  geom_smooth(formula=y~s(x,k=20),method="gam") + theme_classic()+
  scale_fill_continuous(low="white",high="black", trans="log",
                        breaks=c(1,10,100,1000,10000))

ph <- ggplot(plot.df2, aes(par_hd,hyb)) + geom_hex() + facet_grid(.~type,scales="fixed")+
  theme_classic()+
  scale_fill_continuous(low="white",high="black", trans="log",
                        breaks=c(1,10,100,1000,10000))

