#setwd("R")
source("../../chapter2/R/util/rhat.R")

s1 <- readRDS("arrayOutput/chain1_sb.rds")
class(s1) <- "myMcmcObj"
s2 <- readRDS("arrayOutput/chain2_sb.rds")
class(s2) <- "myMcmcObj"
s3 <- readRDS("arrayOutput/chain3_sb.rds")
class(s3) <- "myMcmcObj"
s4 <- readRDS("arrayOutput/chain4_sb.rds")
class(s4) <- "myMcmcObj"

alphadf <- data.frame(alpha=c(s1$samples$alpha,
                              s2$samples$alpha,
                              s3$samples$alpha,
                              s4$samples$alpha))
alphadf$chain <- factor(rep(1:4,each=10000))
alphadf$iter <- rep(1:10000,4)

library(ggplot2)
ggplot(alphadf, aes(x=iter, y=alpha, color=chain,group=chain)) + 
  geom_line()+theme_classic()+
  scale_color_brewer(palette = "RdYlGn")

ggplot(alphadf, aes(x=alpha, group=chain, color=chain, fill=chain)) + geom_density(alpha=.5)

G <- dim(s1$summaries$means_betas)[2]
upper.bound <- data.frame(lb=4*G*exp(-(2^15-1)/alphadf$alpha))
ggplot(upper.bound, aes(x=lb)) + theme_classic() + scale_x_continuous(trans="log")+
  geom_histogram()
