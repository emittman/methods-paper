roc1 <- readRDS("../data/roc-df1.rds")
roc2 <- readRDS("../data/roc-df2.rds")
roc3 <- readRDS("../data/roc-df3.rds")


rr <- dplyr::filter(roc2, threshold==.25, p==2, type=="BNP")

roc_fns <- rbind(roc1,roc2,roc3) %>%
  dlply(.(threshold, p), function(sss){
    dlply(sss, .(type), function(ss){
      dlply(ss, .(sim), function(s){
        return(approxfun(s$FPR, s$TPR, method = "linear", ties = max))
      })
    })
  })

grid <- 1:100*.001

mean.rocs <- ldply(roc_fns, function(setting){
  ldply(setting, function(type){
    all <- ldply(type, function(f_sim){
      data.frame(FPR=grid, TPR=f_sim(grid))
    })
    ddply(all, .(FPR), summarise, 
          lower=quantile(TPR, .25,na.rm=TRUE),
          upper=quantile(TPR, .75, na.rm=TRUE),
          TPR = mean(TPR))
  })
})

str(mean.rocs)

mean.rocs$p <- factor(mean.rocs$p, levels=2:5, labels=c("Parental HD", "Hybrid","Hybrid HD", "Flow Cell"))

ggplot(mean.rocs, aes(FPR, TPR, color=type)) + geom_line() + 
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=type), alpha=.5)+
  facet_grid(threshold~p) + theme_bw(base_size=14)
