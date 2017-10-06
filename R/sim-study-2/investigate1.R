library(limma)

load("R/data/heterosis_design.RData")

library(edgeR)
library(dplyr)
dge <- DGEList(sim$y) %>% calcNormFactors(method="TMM")

voom_dge <- voom(dge, design = X)
fit <- lmFit(voom_dge)
df.voom <- data.frame(coef(fit))
names(df.voom)[1] <- "intercept"

df.bnp <- data.frame(t(s$summaries$means_betas))
names(df.bnp) <- names(df.voom)

df.voom$method <- "voom"
df.bnp$method <- "bnp"


library(ggplot2)
library(tidyr)
rbind(df.voom,df.bnp) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=parent_hd, y=hybrid, color=method)) + geom_density_2d(alpha=.2, n=100, binwidth=.01) +
  facet_wrap(~method) + theme_bw()

rbind(df.voom,df.bnp) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=intercept, y=hybrid, color=method)) + geom_density_2d(alpha=.2, n=100, binwidth=.01) +
  facet_wrap(~method) + theme_bw()

rbind(df.voom,df.bnp) %>%
  # gather(key=par, value=value, -method) %>%
  ggplot(aes(x=parent_hd, y=flow_cell, color=method)) + geom_density_2d(alpha=.5, n=50, binwidth=.04) +
  facet_wrap(~method) + theme_bw()

