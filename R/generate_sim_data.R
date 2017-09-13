#truth based on P
s1 <- readRDS("arrayOutput/chain1_sb.rds")
ggplot(data.frame(s1$samples$P[,,50]), aes(pi)) + geom_histogram() +
  scale_x_continuous(trans="log")

pdraw_P <- data.frame(s1$sample$P[,,50])

ggplot(pdraw_P, aes(x=beta.2.,y=beta.3.)) + geom_hex(aes(weight=pi), bins=30)+
  scale_fill_continuous(trans="log",breaks=10^(-10:-1), low="white",high="black")+
  theme_bw()

str(pdraw_P)

# K <- max(s1$samples$num_occupied)

#alphahat <- mean(s1$samples$alpha)
G <- dim(s1$summaries$means_betas)[2]
K <- nrow(pdraw_P)
#sum(alphahat/(alphahat + 1:G-1))
true.zeta <- sample(K, size= G, replace=T, prob = pdraw_P$pi)

x <- length(unique(true.zeta)) %/% 512

maximum.number.of.clusters <- (x+1)*512

true.params <- pdraw_P[true.zeta,]

true.betas <- t(true.params[,2:6])
true.sigmas <- true.params[,7]
load("R/data/heterosis_long.RData")

library(dplyr)
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

N <- nrow(X)
y <- t(sapply(1:G, function(g) rnorm(N, X%*%true.betas[,g], true.sigmas[g])))

