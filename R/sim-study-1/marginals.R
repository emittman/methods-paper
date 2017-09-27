#truth vs. estimates -- pick a sim

i <- sample(10, 1)
mcmc <- readRDS(paste("saved_mcmc/mcmc_sim_",i,sep=""))
sim <- readRDS(paste("sim_data/sim_",i,sep=""))

betas <- data.frame(rbind(cbind(t(mcmc$summaries$means_betas),mcmc$summaries$means_sigmas),
                          cbind(sim$truth$beta, sim$truth$sigma)), check.rows = FALSE, row.names = NULL)
names(betas) <- c(sapply(1:5, function(j) paste("beta[",j,"]",sep="")),"sigma")
betas$type <- rep(c("est","true"), each=10000)

library(tidyr)

betas_long <- gather(betas, key=par, value=value, -type)

source("../../../../cuda_rpackage/R/data.R")
load("../data/heterosis_long.RData")
library(dplyr)
X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

dat <- formatData(counts=sim$y, X=X, transform_y=identity, voom = FALSE)
ind_est <- indEstimates(dat)
prior <- formatPriors(12, estimates=ind_est)

prior_dens <- rbind(data.frame(par="beta[1]",
                               value=rnorm(10000, prior$mu_0[1], 1/sqrt(prior$lambda2[1]))),
                    data.frame(par="beta[2]",
                               value=rnorm(10000, prior$mu_0[2], 1/sqrt(prior$lambda2[2]))),
                    data.frame(par="beta[3]",
                               value=rnorm(10000, prior$mu_0[3], 1/sqrt(prior$lambda2[3]))),
                    data.frame(par="beta[4]",
                               value=rnorm(10000, prior$mu_0[4], 1/sqrt(prior$lambda2[4]))),
                    data.frame(par="beta[5]",
                               value=rnorm(10000, prior$mu_0[5], 1/sqrt(prior$lambda2[5]))),
                    data.frame(par="sigma",
                               value=rgamma(10000, prior$a, prior$b)))
prior_dens$type <- "G[0]"
library(ggplot2)
p <- ggplot(rbind(betas_long,prior_dens), aes(x=value)) +
  geom_density(aes(group=type, fill=type, linetype=type), alpha=.5) +
  facet_wrap(~par, scales="free", labeller = label_parsed, nrow = 3)+
  # scale_y_continuous(trans="log", breaks=c(10^(-5:2)))+
  scale_fill_manual(values=c("est"="blue","G[0]"="green","true"="red"),
                     labels=c("Estimate", expression(G[0]),"True value"))+
  scale_linetype_manual(values=c("est"=1,"G[0]"=3,"true"=2),
                    labels=c("Estimate", expression(G[0]),"True value"))+
  theme_bw(base_size=14) + xlab("")
  

p
ggsave("../../figures_tables/par-marginal-sim-1.pdf", width=8, height=8)



