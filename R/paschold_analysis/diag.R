setwd("R/arrayOutput/")
source("../../../chapter2/R/util/rhat.R")

s <- readRDS("chain1_sb_1024.rds")
s$samples$beta <- NULL
s1 <- s

s <- readRDS("chain2_sb_1024.rds")
s$samples$beta <- NULL
s2 <- s

s <- readRDS("chain3_sb_1024.rds")
s$samples$beta <- NULL
s3 <- s

s <- readRDS("chain4_sb_1024.rds")
s$samples$beta <- NULL
s4 <- s

gf <- gelman.factors(s1,s2,s3,s4,n_iter = 1000)
str(gf)

library(ggplot2)
rhatdf <- with(gf, data.frame(g=rep(1:36821, 5), p=rep(1:5, each=36821), Rhat = as.numeric(t(beta))))
               
ggplot(rhatdf, aes(x=Rhat)) + facet_wrap(~p, nrow=1, scales="free") + geom_histogram()

summary(gf$sigma)

summaries <- combine_chains_summaries(s1,s2,s3,s4)

library(tidyr)
prob.df <- data.frame(t(summaries$prob))
C <- list(extreme_heterosis = matrix(c(0, 1, 1, 1, 0,
                                       0, 1, 1,-1, 0,
                                       0,-1, 1, 1, 0,
                                       0,-1, 1,-1, 0),4, 5, byrow=T),
          high_mean   = matrix(c(0,1,1,0,0,
                                 0,-1,1,0,0), 2,5,byrow=T),
          hp_h12      = matrix(c(0,1,1,1,0,
                                 0,-1,1,1,0), 2,5,byrow=T),
          lp_h12      = matrix(c(0,1,-1,-1,0,
                                 0,-1,-1,-1,0), 2, 5, byrow=T),
          hp_h21      = matrix(c(0,1,1,-1,0,
                                 0,-1,1,-1,0), 2,5, byrow=T),
          lp_h21      = matrix(c(0,1,-1,1,0,
                                 0,-1,-1,1,0), 2, 5, byrow=T),
          de_p1       = matrix(c(0,1,0,0,0), 1, 5, byrow=T),
          high_or_mid = matrix(c(0,0,1,0,0), 1, 5, byrow=T))
names(prob.df) <- names(C)
gather(prob.df, key=type, value=value) %>%
  ggplot(aes(x=value)) + geom_histogram(bins=100) +
  facet_wrap(~type, scales = "free") + 
  scale_y_continuous(trans="log", breaks=c(1,10,100,1000))


#compare cis based on moment matching vs. full posterior
s <- readRDS("chain1_sb_1024.rds")
id <- sort(sample(1:36000, 5))
samp_samp <- s$samples$beta[,,id]
rm(list="s")
library(plyr)
cis <- adply(samp_samp, c(1,3), function(x){
  q <- quantile(x, c(.025,.975))
  data.frame(lower=q[1],upper=q[2],type="full bayes")
})

cis_asymp <- ldply(id, function(g){
  mu <- summaries$means_betas[,g]
  es <- summaries$meansquares_betas[,g]
  sigma <- sqrt(es - mu^2)
  data.frame(v=1:5, g=g, lower=mu-1.96*sigma, upper=mu+1.96*sigma, type="asymptotic")
})


rbind(cis, cis_asymp) %>%
  mutate(type = factor(type)) %>%
  ggplot(aes(x=type,xend=type, y=lower,yend=upper, lty=type))+
  geom_segment() + facet_wrap(g~v, scales="free")+coord_flip()
