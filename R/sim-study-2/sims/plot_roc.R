roc1 <- readRDS("../data/roc-df1.rds")
roc2 <- readRDS("../data/roc-df2.rds")
roc3 <- readRDS("../data/roc-df3.rds")

consolidated <- rbind(roc1,roc2,roc3) %>%
  group_by(sim, threshold, p, type, FPR) %>%
  summarize(TPR = max(TPR)) %>%
  mutate(id = paste(type,"_",sim,sep=""))

consolidated$p <- factor(consolidated$p, levels=2:5,
                         labels=c("parental HD", "hybrid",
                                  "hybrid HD", "flow cell"))
ssss <- sample(1:10,1)
p1 <- filter(consolidated,threshold<.8, threshold>0, sim==ssss, type != "voom-limma") %>%#,p %in% c("parental HD","hybrid","hybrid HD")) %>%
  # ddply(.(threshold,p,type,FPR), summarise, TPR=mean(TPR)) %>%
  ggplot(aes(FPR,TPR, color=type, linetype=type, group=id)) + geom_line() + 
  facet_grid(threshold~p, scales = "free")+theme_bw(base_size=14)+
  theme(legend.position = c(.9,.14),
        legend.margin = margin(-10,-10,-10,-10,unit="pt"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))
p1 <- p1 + theme(axis.text.x = element_text(angle = 90, vjust =.5))
aucs <- arrange(consolidated, sim, threshold, p, type, FPR) %>%
  filter(FPR < .1, type!="voom-limma") %>%
  ddply(.(sim, threshold, p, type), function(x){
    width = diff(x$FPR)
    height = (x$TPR[-1] + x$TPR[-length(x$TPR)])/2
    data.frame(AUC = sum(width*height))
  })

p2 <- filter(aucs, threshold>0,threshold<.8)%>%#, p %in% c("parental HD","hybrid","hybrid HD")) %>%
  ggplot(aes(x=type, y=AUC)) + geom_boxplot(aes(color=type)) +
  geom_line(aes(group=sim), alpha=.5) + facet_grid(threshold~p) +
  theme_bw(base_size=14)+xlab("")+
  theme(legend.position = "none",
  axis.text.x = element_text(angle = 90, vjust = .5))

library(cowplot)
plot_grid(p1,p2,ncol=1)
ggsave("../../../figures_tables/ss2-roc-auc.pdf", width=6, height=10)

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
mean.rocs %>% filter(type != "voom-limma") %>%
ggplot(aes(FPR, TPR, color=type)) + geom_line() + 
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=type), alpha=.5)+
  facet_grid(threshold~p) + theme_bw(base_size=14) +
  scale_x_continuous(limits=(c(.001, .05)), breaks=c(0.0,.02,.04))
