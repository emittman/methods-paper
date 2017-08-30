library(Paschold2012)
library(dplyr)
library(tidyr)
# source("getNormFactors.R")

#differential expression, B73 v. B73xMo17
my_dat_wide <- Paschold2012 %>%
  mutate(genotype_replicate = paste(genotype,replicate,sep="_")) %>%
  select(GeneID, genotype_replicate, total) %>%
  tidyr::spread(genotype_replicate, total) %>%
  select(1:17)

#remove all zero genes
zero_id <- which(apply(data.matrix(my_dat_wide[,-1]), 1, function(x) all(x==0)))
my_dat_wide <- my_dat_wide[-zero_id,]

#transform counts to log (+1) scale, normalize
# my_dat_wide[,2:17] <- normalizeData(data.matrix(my_dat_wide[,2:17]),
                                   # group=rep(1:4, each=4), trans.to.log = TRUE)

save(my_dat_wide, file="data/heterosis_counts.RData")

#convert to long format, add indicators for genotype and flowcell
my_dat <- tidyr::gather(my_dat_wide, key=genotype_replicate, value=total, -GeneID) %>%
  extract(col=genotype_replicate, into=c("genotype","replicate"),
                     regex = "([[:alnum:]]+)_([[:alnum:]]+)", convert=TRUE) %>%
  mutate(parent_hd = ifelse(genotype=="B73",1, ifelse(genotype=="Mo17",-1,0)),
         hybrid = ifelse(genotype %in% c("B73xMo17","Mo17xB73"),1,0),
         hybrid_hd = ifelse(genotype %in% c("B73","Mo17"), 0, ifelse(genotype == "B73xMo17", 1, -1)),
         flow_cell = ifelse(replicate %in% c(1,2), 1, -1)) %>%
  arrange(GeneID)

save(my_dat, file="data/heterosis_long.RData")

source("../../../cuda_rpackage/R/data.R")

X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>%
  model.matrix(~parent_hd+hybrid+hybrid_hd+flow_cell, data=.)

y <- data.matrix(my_dat_wide[,2:17])

cuda_dat <- formatData(counts = y, groups = rep(1:2, each=4), X = X, voom = TRUE, transform_y=identity)
ind_est <- indEstimates(cuda_dat)

save(ind_est, file="data/ind-est-heterosis.RData")
save(cuda_dat, file="data/cuda-dat-heterosis.Rdata")


# ##voom fiddle
# library(limma)
# obj <- voomWithQualityWeights(y, X, plot = TRUE)
# fit <- lmFit(object = obj)
# library(ggplot2)
# library(GGally)
# my_hex <- function(data, mapping, ..., low = "white", high ="blueviolet") {
#   ggplot(data = data, mapping = mapping) +
#     geom_hex(...) +
#     scale_fill_gradient(low = low, high = high, trans="log", breaks=c(1,10,100,1000,10000))+theme_bw()
# }
# ggpairs(data.frame(coefficients(fit)), lower = list(continuous=my_hex))
# 
# voom_dat <- formatData(my_dat_wide[,2:17], X = X, voom=TRUE)
# voom_est <- indEstimates(voom_dat)
# data.frame(intercept=voom_est$beta[1,], log_sigma=.5*log(voom_est$sigma2)) %>%
#   ggplot(aes(x=intercept,y=log_sigma)) + geom_hex()     
#   scale_fill_gradient(low = "white", high = "blueviolet", trans="log", breaks=c(1,10,100,1000,10000))+
#   theme_bw()
# 
# plot(voom_est$beta[1,], sqrt(voom_est$sigma2))
