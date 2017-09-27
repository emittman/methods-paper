#geweke

#sim1
mcmc1 <- readRDS("mcmc_sim_1")
str(mcmc1$samples$beta)

library(plyr)
library(coda)

geweke <- adply(mcmc1$samples$beta, c(1, 3), function(ch){
  data.frame(geweke = geweke.diag(ch)[[1]])
})

library(ggplot2)

ggplot(geweke, aes(x=geweke)) + geom_density(fill="blue", alpha=.5) +
  facet_wrap(~v) + stat_function(fun = dnorm)
