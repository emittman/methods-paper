set.seed(101117345)

load("R/data/heterosis_design.RData")

dir1 <- sample(pi/4 + c(0, pi/2), 10000, replace=TRUE)
dir2 <- sample(pi/4 + c(0, pi/2), 10000, replace=TRUE)
dir3 <- sample(pi/4 + c(0, pi/2), 10000, replace=TRUE)
dir4 <- sample(pi/4 + c(0, pi/2), 10000, replace=TRUE)
dir5 <- sample(pi/4 + c(0, pi/2), 10000, replace=TRUE)

mag <- rnorm(10000, 2, .5)

x1 <- mag*cos(dir1) + 2.5
x2 <- mag*sin(dir1)*cos(dir2)
x3 <- mag*sin(dir1)*sin(dir2)*cos(dir3)
x4 <- mag*sin(dir1)*sin(dir2)*sin(dir3)*cos(dir4)
x5 <- mag*sin(dir1)*sin(dir2)*sin(dir3)*sin(dir4)*cos(dir5)
x6 <- (mag*sin(dir1)*sin(dir2)*sin(dir3)*sin(dir4)*sin(dir5) - .8)

plot(x2,x4)
abline(0,1)
abline(0,-1)

df <- data.frame(x1,x2,x3,x4,x5,exp(x6))

X
dat <- t(apply(data.matrix(df), 1, function(x){
  rnorm(16, X%*%x[1:5], x[6])
}))

head(dat)
est <- t(apply(dat, 1, function(y){
  c(coef(lm(y~0+X)), summary(lm(y~0+X))$sigma)
}))

pairs(df[,1:5], pch='.')
library(ggplot2)
library(GGally)
my_hex <- function(data, mapping, ...) {
    ggplot(data = data, mapping=mapping) +
    geom_hex(aes(fill=..density..), ..., bins=30, color = NA)+
    scale_fill_continuous(trans="log",low="white",high="darkviolet")+
    geom_abline(slope=c(-1,1))
}
p <- ggpairs(df, lower = list(continuous=my_hex))
q <- ggpairs(as.data.frame(est), lower= list(continuous=my_hex))

y <- as.data.frame(round(exp(dat)))
names(y) <- c('a_1','a_2','a_3','a_4',
              'ab_1','ab_2','ab_3','ab_4',
              'b_1','b_2','b_3','b_4',
              'ba_1','ba_2','ba_3','ba_4')
head(y)

setwd("R/sim-study-3/dat")
save(y, file="R/sim-study-3/data/counts_1.RData")
save(X, file="R/sim-study-3/data/design.RData")
save(df, file="R/sim-study-3/data/truth.RData")
# library(tidyr)
# library(dplyr)
# df_y <- mutate(y, gene = 1:10000)
# df_y <- df_y %>%
#   gather(key=key,value=total_count,1:16) %>%
#   extract(col = key, into = c("geno","rep"), remove = TRUE,
#           regex = "^([[:alnum:]]+)_([[:alnum:]]+)$")
# 
# filter(df_y, gene<5) %>%
#   ggplot(aes(x=geno, y=total_count+1)) + geom_boxplot() + facet_wrap(~gene, scales="free", ncol=1) +
#   scale_y_continuous(trans="log")
