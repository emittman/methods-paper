#reduce beta draws... large files may require a restart of R

s <- readRDS("../arrayOutput/chain1_sb_1024.rds")
with(s$samples, saveRDS(list(beta[,1:250*4,],
                             tau2[,1:250*4],
                             alpha),
                        file="chain1samples.rds"))
s <- readRDS("../arrayOutput/chain2_sb_1024.rds")
with(s$samples, saveRDS(list(beta[,1:250*4,],
                             tau2[,1:250*4],
                             alpha),
                        file="chain2samples.rds"))
s <- readRDS("../arrayOutput/chain3_sb_1024.rds")
with(s$samples, saveRDS(list(beta[,1:250*4,],
                             tau2[,1:250*4],
                             alpha),
                        file="chain3samples.rds"))
s <- readRDS("../arrayOutput/chain4_sb_1024.rds")
with(s$samples, saveRDS(list(beta[,1:250*4,],
                             tau2[,1:250*4],
                             alpha),
                        file="chain4samples.rds"))
#clear workspace
library(abind)
s <- readRDS("chain1samples.rds")
s2 <- readRDS("chain2samples.rds")
names(s) <- c("beta","tau2","alpha")
s$beta <- abind(s$beta, s2[[1]], along = 2)
s$tau2 <- abind(s$tau2, s2[[2]], along = 2)
s$alpha <- c(s$alpha, s2[[3]])
rm(list="s2")

s3 <- readRDS("chain3samples.rds")
s$beta <- abind(s$beta, s3[[1]], along = 2)
s$tau2 <- abind(s$tau2, s3[[2]], along = 2)
s$alpha <- c(s$alpha, s3[[3]])
rm(list="s3")

s4 <- readRDS("chain4samples.rds")
s$beta <- abind(s$beta, s4[[1]], along = 2)
s$tau2 <- abind(s$tau2, s4[[2]], along = 2)
s$alpha <- c(s$alpha, s4[[3]])
rm(list="s4")

saveRDS(s, file="samples.rds")

load("../data/heterosis_design.RData")
mu_B <- sapply(1:36821, function(g) X[1,] %*% s$beta[,,g])
mu_BM <- sapply(1:36821, function(g) X[5,] %*% s$beta[,,g])
mu_M <- sapply(1:36821, function(g) X[9,] %*% s$beta[,,g])
mu_MB <- sapply(1:36821, function(g) X[13,] %*% s$beta[,,g])

b <- mean(mu_B[,id])
bm <- mean(mu_BM[,id])
m <- mean(mu_M[,id])
mb <- mean(mu_MB[,id])

bq <- quantile(mu_B[,id], c(.05,.95))
bmq <- quantile(mu_BM[,id], c(.05,.95))
mq <- quantile(mu_M[,id], c(.05,.95))
mbq <- quantile(mu_MB[,id], c(.05,.95))

#From compare_orig.R
sv2+geom_pointrange(data=data.frame(genotype=c("B73","B73xMo17","Mo17xB73","Mo17"),
                                    total = c(b, bm, mb, m), lower=c(bq[1],bmq[1],mbq[1],mq[1]),
                                    upper=c(bq[2],bmq[2],mbq[2],mq[2])), aes(ymin=lower,ymax=upper), color="red")


th <- .05
high_parent_th12 <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,1,0)%*%s$beta[,,g] > th),
    as.numeric(c(0, 1,1,1,0)%*%s$beta[,,g] > th))
}), 2, mean)

high_parent_th21 <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,-1,0)%*%s$beta[,,g] > th),
       as.numeric(c(0, 1,1,-1,0)%*%s$beta[,,g] > th))
}), 2, mean)

high_parent_mean <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,0,0)%*%s$beta[,,g] > th),
       as.numeric(c(0, 1,1,0,0)%*%s$beta[,,g] > th))
}), 2, mean)

low_parent_mean <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,0,0)%*%s$beta[,,g] < -th),
       as.numeric(c(0, 1,1,0,0)%*%s$beta[,,g] < -th))
}), 2, mean)

low_parent_th12 <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,1,0)%*%s$beta[,,g] < -th),
       as.numeric(c(0, 1,1,1,0)%*%s$beta[,,g] < -th))
}), 2, mean)

lows_parent_th21 <- apply(sapply(1:36821, function(g){
  pmin(as.numeric(c(0,-1,1,-1,0)%*%s$beta[,,g] < -th),
       as.numeric(c(0, 1,1,-1,0)%*%s$beta[,,g] < -th))
}), 2, mean)

high_mean_effsize <- apply(sapply(1:36821, function(g){
  pmax(pmin(as.numeric(c(0,-1,1,0,0)%*%s$beta[,,g]),
       as.numeric(c(0, 1,1,0,0)%*%s$beta[,,g])),0)
}), 2, mean)
