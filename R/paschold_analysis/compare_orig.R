#discovery vs nondiscovery

library(readxl)

pasch <- read_excel("paschold_analysis/Copy of TableS3.xlsx", sheet = 1, skip = 1)

#remove rows containing summaries
pasch <- pasch[1:39656,]

table(pasch$`MxB class`)
table(pasch$`BxM class`)

#\citet{paschold} determined a classification of genes by combining information from
# several hypothesis tests. In Figure \ref{orig-compare}, we show histograms of posterior 
#probabilities obtained by our method split by the classes as reported in \citet{paschold}

hp21.GeneID <- pasch$Gene[which(pasch$`MxB class` %in% c(5,6))]
lp21.GeneID <- pasch$Gene[which(pasch$`MxB class` %in% c(7,8))]
hp12.GeneID <- pasch$Gene[which(pasch$`BxM class` %in% c(5,6))]
lp12.GeneID <- pasch$Gene[which(pasch$`BxM class` %in% c(7,8))]

#load subsetted data with GeneID
load("data/heterosis_counts.RData")

hp21.id <- match(hp21.GeneID, my_dat_wide$GeneID)
lp21.id <- match(lp21.GeneID, my_dat_wide$GeneID)
hp12.id <- match(hp12.GeneID, my_dat_wide$GeneID)
lp12.id <- match(lp12.GeneID, my_dat_wide$GeneID)

load("data/heterosis_design.RData")
load("paschold_analysis/contrasts.RData")
load("paschold_analysis/combined_summaries.RData")
probs.df <- data.frame(t(summaries$prob))
names(probs.df) <- names(C)

probs.df$hp12 <- FALSE
probs.df$hp12[hp12.id] <- TRUE
probs.df$hp21 <- FALSE
probs.df$hp21[hp21.id] <- TRUE
probs.df$lp12 <- FALSE
probs.df$lp12[lp12.id] <- TRUE
probs.df$lp21 <- FALSE
probs.df$lp21[lp21.id] <- TRUE

#Using a threshold (see combine_samples.R)
probs.df$high_mean <- high_parent_mean
probs.df$hp_h12 <- high_parent_th12
probs.df$lp_h12 <- low_parent_th12
probs.df$hp_h21 <- high_parent_th21
probs.df$lp_h21 <- low_parent_th21

library(ggplot2)
library(tidyr)
library(dplyr)

plot.df <- select(probs.df, hp_h12, hp_h21, lp_h12, lp_h21, hp12, hp21, lp12, lp21) %>%
  gather(key=`hypothesis`, value=`posterior probability`, hp_h12, lp_h21, lp_h12, hp_h21) %>%
  mutate(discovery = factor(ifelse(hypothesis=='hp_h12', hp12,
                                   ifelse(hypothesis=='hp_h21', hp21,
                                          ifelse(hypothesis=='lp_h12', lp12, lp21))),
                            levels=c(T,F),
                            labels=c("Discovery","Non-discovery"))) %>%
  mutate(hypothesis = factor(hypothesis,
                             levels = c("hp_h12","hp_h21","lp_h12","lp_h21"),
                             labels = c("high-parent B73xMo17",
                                        "high-parent Mo17xB73",
                                        "low-parent B73xMo17",
                                        "low-parent Mo17xB73")))

ggplot(plot.df, aes(x=`posterior probability`)) +
  geom_histogram(bins=30) + facet_grid(discovery ~ hypothesis, scales="free") +
  theme_bw()

#Will's results lots of cross-contamination from combind_samples.R
wills <- read.csv("paschold_analysis/data.csv")
wills <- filter(wills, geneID %in% my_dat_wide$GeneID) %>%
  select(geneID, c(18,20:23))
names(wills)[-1] <- names(probs.df)[2:6]
wills$type <- "parametric"
mine <- probs.df[,2:6]
mine$geneID <- wills$geneID
mine$type <- "BNP"
compare.df <- rbind(wills, mine) %>% 
  gather(key=`heterosis type`, value=value, high_mean, hp_h12, lp_h12, hp_h21, lp_h21) %>%
  spread(key=type, value=value)
  
ggplot(data=compare.df, aes(x=parametric, y=BNP)) + geom_hex(bins=30) + facet_wrap(~`heterosis type`) +
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+theme_bw()+
  geom_abline(slope=1) + coord_fixed() + geom_smooth()

bnp_high.param_low <- compare.df$geneID[order(compare.df$BNP-compare.df$parametric,decreasing=T)[1:3]]
my_dat_wide[my_dat_wide$GeneID %in% bnp_high.param_low, ]
compare.df[which(compare.df$geneID %in% bnp_high.param_low),]

load("data/heterosis_long.RData")
load("data/cuda-dat-heterosis.Rdata")
#normalized log-counts in long format
my_dat2 <- data.frame(data.matrix(cuda_dat$transformed_counts)) %>%
  mutate(GeneID= my_dat_wide$GeneID) %>%
  gather(key=geno_rep, value= total, -GeneID) %>%
  extract(geno_rep, into = c("genotype", "replicate"),
          regex = "^([[:alnum:]]+)_([[:alnum:]]+)$",
          convert = TRUE) %>%
  mutate(genotype = factor(genotype,
                levels=c("B73","B73xMo17","Mo17xB73","Mo17")))

sv1 <- filter(my_dat2, GeneID %in% bnp_high.param_low[1]) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1) +
  facet_wrap(~GeneID)

bnp_low.param_high <- compare.df$geneID[order(-compare.df$BNP+compare.df$parametric, decreasing=TRUE)[1:3]]
my_dat_wide[my_dat_wide$GeneID %in% bnp_low.param_high, ]
compare.df[which(compare.df$geneID %in% bnp_low.param_high),]


sv2 <- filter(my_dat2, GeneID %in% bnp_low.param_high[1]) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1)# + facet_wrap(~GeneID)+
  scale_y_continuous(trans="log", breaks=2^1:15)


cowplot::plot_grid(sv1+ggtitle("Gene X"),
                   sv2+ggtitle("Gene Y"), nrow=2)

both.high <- compare.df$geneID[order(compare.df$BNP+compare.df$parametric,decreasing=T)[1:3]]
my_dat_wide[my_dat_wide$GeneID %in% both.high, ]
compare.df[which(compare.df$geneID %in% both.high),]
filter(my_dat, GeneID %in% both.high) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1) +
  facet_wrap(~GeneID, scales="free") +
  scale_y_continuous(trans="log")

coin.flip <- compare.df$geneID[order(abs(compare.df$BNP-1)+abs(compare.df$parametric-1), decreasing=FALSE)[1:10]]
compare.df[which(compare.df$geneID %in% coin.flip),]
r3 <- my_dat_wide[my_dat_wide$GeneID %in% coin.flip, ]
filter(my_dat, GeneID %in% coin.flip[1]) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1,height=0)
  # facet_wrap(~GeneID, scales="free") +
  # scale_y_continuous(trans="log")


effsize_par <- read.csv("paschold_analysis/data.csv")[,c(1,36)]
effsize_par <- filter(effsize_par, geneID %in% my_dat_wide$GeneID)
effsize_par$BNP <- high_mean_effsize
ggplot(effsize_par,
       aes(high.parent_hybrid.mean_effect_size, BNP)) +
  geom_hex(aes(fill=..density..), bins=30)+
  scale_fill_continuous(trans="log",low="white",high="midnightblue")+
  theme_bw()

both.high <- filter(compare.df, `heterosis type` == "high_mean") %>%
  mutate(comb.rank = BNP + parametric) %>%
  arrange(-comb.rank)

both.low <- filter(compare.df, `heterosis type` %in% c("lp_h12","lp_h21")) %>%
  group_by(geneID) %>%
  summarise(comb.rank=sum(BNP + parametric)) %>%
  arrange(-comb.rank)

opposite <- filter(compare.df, `heterosis type` %in% c("lp_h12","hp_h21")) %>%
  group_by(geneID) %>%
  summarise(comb.rank=sum(BNP + parametric)) %>%
  arrange(-comb.rank)

mids <- compare.df %>%
  group_by(geneID) %>%
  summarise(comb.rank=sum(BNP + parametric)) %>%
  arrange(comb.rank)

#exemplar gene ids
ids <- c(as.character(both.high$geneID[1]),as.character(both.low$geneID[1]),
         as.character(opposite$geneID[1]), as.character(mids$geneID[1]))

idx <- match(ids, my_dat_wide$GeneID)

s <- readRDS("paschold_analysis/samples.rds")
mu_B <- sapply(1:36821, function(g) X[1,] %*% s$beta[,,g])
mu_BM <- sapply(1:36821, function(g) X[5,] %*% s$beta[,,g])
mu_M <- sapply(1:36821, function(g) X[9,] %*% s$beta[,,g])
mu_MB <- sapply(1:36821, function(g) X[13,] %*% s$beta[,,g])

cis <- plyr::ldply(idx, function(id){
  b <- mean(mu_B[,id])
  bm <- mean(mu_BM[,id])
  m <- mean(mu_M[,id])
  mb <- mean(mu_MB[,id])
  
  bq <- quantile(mu_B[,id], c(.05,.95))
  bmq <- quantile(mu_BM[,id], c(.05,.95))
  mq <- quantile(mu_M[,id], c(.05,.95))
  mbq <- quantile(mu_MB[,id], c(.05,.95))
  data.frame(geneID=my_dat_wide$GeneID[id],
             genotype=c("B73","B73xMo17","Mo17xB73","Mo17"),
             total = c(b, bm, mb, m),
             lower=c(bq[1],bmq[1],mbq[1],mq[1]),
             upper=c(bq[2],bmq[2],mbq[2],mq[2]))
})
cis <- mutate(cis,
              geneID = factor(geneID, levels=ids))#, labels = 
                                #c("HPH", "LPH", "mixed", "neither")))


filter(my_dat2, GeneID %in% ids) %>%
  mutate(geneID = factor(GeneID, levels=ids))%>%#, labels = 
                           # c("HPH", "LPH", "mixed", "neither"))) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1) +
  geom_pointrange(data=cis, aes(ymin=lower,ymax=upper), color="red") +
  facet_wrap(~geneID, scales = "free_y") + theme_bw(base_size=14) + 
  xlab("") + ylab("") +
  ggtitle("Estimated mean log cpm by genotype")+
  theme(axis.text.x = element_text(angle = 15, hjust = .9, vjust=.9))

library(xtable)
select_genes <- my_dat_wide[idx,]
select_genes$`example` <- c("HPH","LPH", "mixed","neither")
str(select_genes)
select(select_genes, 18,2:9,14:17,10:13) %>%
  xtable %>%
  print.xtable(include.rownames = FALSE, booktabs=TRUE, scalebox=.9)

# 2-way table of binned posterior probs
new <- compare.df
new$dBNP<- cut(compare.df$BNP, breaks = c(0,.05,.15,.85,.95,1),
               include.lowest = TRUE, right=TRUE)
new$dparametric <- cut(compare.df$parametric, breaks = c(0,.05,.15,.85,.95,1),
                       include.lowest = TRUE, right=TRUE)

make_margin_table <- function(table){
  cbind(rbind(table/margin.table(table),
              total=margin.table(table,2)/margin.table(table)),
        total=c(margin.table(table,1)/margin.table(table),1))
}

df_lph12 <- filter(new, `heterosis type`=="lp_h12")
tbl_lph12 <-   table(df_lph12$dBNP, df_lph12$dparametric)
print.xtable(xtable(make_margin_table(tbl_lph12), digits=3),
             include.rownames = T,booktabs = TRUE)
df_lph21 <- filter(new, `heterosis type`=="lp_h21")
tbl_lph21 <-   table(df_lph21$dBNP, df_lph21$dparametric)
print.xtable(xtable(make_margin_table(tbl_lph21), digits=3),
             include.rownames = T,booktabs = TRUE)

df_hph12 <- filter(new, `heterosis type`=="hp_h12")
tbl_hph12 <-   table(df_hph12$dBNP, df_hph12$dparametric)
print.xtable(xtable(make_margin_table(tbl_hph12), digits=3),
             include.rownames = T,booktabs = TRUE)

df_hph21 <- filter(new, `heterosis type`=="hp_h21")
tbl_hph21 <-   table(df_hph21$dBNP, df_hph21$dparametric)
print.xtable(xtable(make_margin_table(tbl_hph21), digits=3),
             include.rownames = T,booktabs = TRUE)

ggplot(data=new, aes(x=dparametric, y=dBNP)) + geom_hex(bins=30) + facet_wrap(~`heterosis type`) +
scale_fill_continuous(trans="log",low="white",high="midnightblue")+theme_bw()+
geom_abline(slope=1) + coord_fixed() + geom_smooth()

biggerthan3 <- my_dat_wide$GeneID[which(summaries$means_betas[1,]>3)]
filter(new) %>%
  with(data=., table(dBNP,dparametric,`heterosis type`)) #%>%
  xtable() %>%
  print()


oe_id <- new$geneID#[which(new$dBNP == "(0,0.25]" & new$dparametric == "(0.75,1]")]
set.seed(12121981)
oe_ids <- match(sample(oe_id, 4),my_dat_wide$GeneID)

cis_beta_full <- plyr::ldply(oe_ids, function(g){
  q <- apply(s$beta[,,g], 1, function(p) quantile(p, c(.05,.95)))
  m <- apply(s$beta[,,g], 1, mean)
  out <- data.frame(g=my_dat_wide$GeneID[g],mean = m, lower = q[1,], upper = q[2,], p=1:5)
})

samp_full <- plyr::ldply(oe_ids, function(g){
  data.frame(g=my_dat_wide$GeneID[g],value=c(t(s$beta[,,g])), p=rep(1:5, each=1000))
})

cis_beta_approx <- plyr::ldply(oe_ids, function(g){
  m = summaries$means_betas[,g]
  s = with(summaries, sqrt(meansquares_betas[,g] - means_betas[,g]^2))
  lower = m - qnorm(.95)*s
  upper = m + qnorm(.95)*s
  data.frame(g=my_dat_wide$GeneID[g], mean=m, sd=s, p=1:5, lower=lower, upper=upper)
})

curve_approx <- plyr::ddply(cis_beta_approx, plyr::.(g, p), function(i){
  q=quantile(filter(samp_full, g==i$g, p==i$p)$value, c(0,1))
  x=seq(from=q[1],to=q[2], length.out=100)
  y=dnorm(x, i$mean, i$sd)
  data.frame(x=x, value=y)
})

yscales <- plyr::ddply(samp_full, plyr::.(g, p), function(i){
  data.frame(lim = max(density(i$value)$y)/6)
})

cis_beta_approx <- mutate(cis_beta_approx,
                          lim = -1*yscales$lim[which(g==g, p==p)])

cis_beta_full <- mutate(cis_beta_full,
                        lim = -.5*yscales$lim[which(g==g, p==p)])

samp_ci_plot <- filter(samp_full) %>%
  mutate(g=factor(g, levels=levels(g)[order(filter(cis_beta_approx, p==1)$mean)])) %>%
ggplot(aes(x=value)) + geom_density(color="red") +
  facet_wrap(g~p, scales="free",nrow=4) +
  geom_segment(data=cis_beta_full, inherit.aes = FALSE,
               aes(x=lower,xend=upper, y=lim, yend=lim), lty="solid", color="red") +
  geom_segment(data=cis_beta_approx, inherit.aes = FALSE,
               aes(x=lower,xend=upper, y=lim, yend=lim), lty="solid", color="black")+
  geom_line(data=curve_approx, inherit.aes = FALSE, aes(x=x, y=value),
            color="black")+ggtitle(expression(paste("Posteriors and 90% CIS for ",beta, ", and normal approximations")))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 25, hjust=.9, vjust=.9))


samp_traceplot <- filter(samp_full)%>%
  mutate(g=factor(g, levels=levels(g)[order(filter(cis_beta_approx, p==1)$mean)])) %>%
  #, g==my_dat_wide$GeneID[oe_ids[4]]) %>%
  mutate(iteration = rep(1:250, 5*4*4),
         chain = factor(rep(rep(1:4, each=250), 5*4))) %>%
  ggplot(aes(x=iteration, y=value, group=chain, color=chain)) + geom_line(alpha=.5) + 
  facet_grid(g~p, scales="free_x")+
  theme_bw() + coord_flip()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust=.9, vjust=.9))+
  geom_hline(data=data.frame(g=rep(unique(samp_full$g),4),
                        p=rep(2:5, each=4),
                        yint=0), aes(yintercept = yint))

my_dnorm <- function(value, type="local") dnorm(value, mean=dfone$mean, sd=dfone$sd)
dfone <- with(cis_beta_approx, data.frame(p=1,g=g[6],mean=mean[6], sd=sd[6]))
dfone2 <- with(dfone, data.frame(value=quantile(filter(samp_full, p==dfone$p, g==dfone$g)$value, c(0,1))))
cbind(dfone,dfone2)

means_cis <- plyr::ldply(oe_ids, function(id){
  b <- mean(mu_B[,id])
  bm <- mean(mu_BM[,id])
  m <- mean(mu_M[,id])
  mb <- mean(mu_MB[,id])
  
  bq <- quantile(mu_B[,id], c(.05,.95))
  bmq <- quantile(mu_BM[,id], c(.05,.95))
  mq <- quantile(mu_M[,id], c(.05,.95))
  mbq <- quantile(mu_MB[,id], c(.05,.95))
  data.frame(geneID=my_dat_wide$GeneID[id],
             genotype=c("B73","B73xMo17","Mo17xB73","Mo17"),
             total = c(b, bm, mb, m),
             lower=c(bq[1],bmq[1],mbq[1],mq[1]),
             upper=c(bq[2],bmq[2],mbq[2],mq[2]))
})
means_cis <- mutate(means_cis,
              geneID = factor(geneID, levels=levels(geneID)[order(filter(cis_beta_approx, p==1)$mean)]))


samp_meanplot <- filter(my_dat2, GeneID %in% my_dat_wide$GeneID[oe_ids]) %>%
  mutate(geneID = factor(GeneID, 
                         levels=levels(cis_beta_approx$g)[order(filter(cis_beta_approx, p==1)$mean)]))%>%#, labels = 
  # c("HPH", "LPH", "mixed", "neither"))) %>%
  ggplot(aes(x=genotype, y=total)) + geom_jitter(width=.1) +
  geom_pointrange(data=means_cis, aes(ymin=lower,ymax=upper), color="red") +
  facet_wrap(~geneID, scales = "free_y") + theme_bw(base_size=14) + 
  xlab("") + ylab("") +
  ggtitle("Estimated mean log cpm by genotype")+
  theme(axis.text.x = element_text(angle = 25, hjust = .9, vjust=.9))
