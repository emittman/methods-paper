#geweke

# #sim1
# mcmc1 <- readRDS("mcmc_sim_1")
# str(mcmc1$samples$beta)

library(plyr)
library(coda)

diagnostics <- function(sims){
  ldply(sims, function(i){
    mcmc <- readRDS(paste("mcmc_sim_",i,sep=""))
    geweke <- adply(mcmc$samples$beta, c(1, 3), function(ch){
      data.frame(geweke = geweke.diag(ch)[[1]],
                 n_eff = effectiveSize(ch))
    })
    geweke
  })
}

diag13 <- diagnostics(1:3)
diag13$sim <- rep(1:3, each=50000)
save(diag13, file="convergence1_3.RData")

diag46 <- diagnostics(4:6)
diag46$sim <- rep(4:6, each=50000)
save(diag46, file="convergence4_6.RData")

diag70 <- diagnostics(7:10)
diag70$sim <- rep(7:10, each=50000)
save(diag70, file="convergence7_10.RData")

library(ggplot2)

rbind(diag13,diag46,diag70) %>%
  mutate(v=factor(v, levels=1:5, labels=c("intercept", "parental HD", "hybrid", "hybrid HD", "flow cell"))) %>%
  ggplot(aes(sample=geweke, group=sim)) + geom_qq(alpha=.5, pch=".") +
  facet_grid(~v) + geom_abline(slope=1, lty=2, color="red")+
  theme_bw()#stat_function(fun = dnorm)

library(dplyr)
mutate(diag13, v=factor(v, levels=1:5, labels=c("intercept", "parental HD", "hybrid", "hybrid HD", "flow cell"))) %>%
ggplot(aes(x=n_eff)) + geom_histogram() +
  facet_grid(sim~v) + scale_x_continuous(breaks=c(0,2000,4000,6000), limits=c(0,6000))

ddply(diag13, .(v), summarise, min_neff = min(n_eff))
