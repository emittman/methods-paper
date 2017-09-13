tanh_trans <- function(scale){
  trans <- function(x) tanh
  inv <- atanh
  scales::trans_new("tanh", transform = trans, inverse = inv)
}

diffHexHist <- function(df1, df2, bins, zmax=NA, zname="A-B"){
  require(hexbin)
  require(ggplot2)
  require(scales)
  stopifnot(dim(df1)[2]==dim(df2)[2]&&dim(df1)[2]==2)
  orig.names <- names(df1)
  names(df1) <- c("x","y")
  names(df2) <- c("x","y")
  makeHexData <- function(df) {
    h <- hexbin(df$x, df$y, xbins=bins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
    data.frame(hcell2xy(h),
               z = h@count/sum(h@count),
               cid = h@cell)
  }  ## find the bounds for the complete data 
  xbnds <- range(c(df1$x, df2$x))
  ybnds <- range(c(df1$y, df2$y))
  
  Ahex <- makeHexData(df1)
  Bhex <- makeHexData(df2)
  
  ##  not all cells are present in each binning, we need to merge by cellID
  byCell <- merge(Ahex, Bhex, by = "cid", all = T)
  
  ##  when calculating the difference empty cells should count as 0
  byCell$z.x[is.na(byCell$z.x)] <- 0
  byCell$z.y[is.na(byCell$z.y)] <- 0
  
  ##  make a "difference" data.frame
  Diff <- data.frame(x = ifelse(is.na(byCell$x.x), byCell$x.y, byCell$x.x),
                     y = ifelse(is.na(byCell$y.x), byCell$y.y, byCell$y.x),
                     z = byCell$z.x - byCell$z.y)
  
  ##  plot the results
  if(is.na(zmax)){
    zmax <- max(abs(Diff))
  }
  
  ggplot(Diff) +
    geom_hex(aes(x = x, y = y, fill = z),
             stat = "identity", alpha = 0.8) +
    scale_fill_gradient2(name=zname, low="blue",high="red", mid="white",
                         trans=scales::trans_new("tanh",
                                                transform = function(x) tanh(x/zmax),
                                                inverse = function(x) atanh(x)*zmax),
                         breaks=c(-zmax, -2:2 * zmax/2)#,
                         #limits=c(-3*zmax,3*zmax)
                         ) +
    scale_x_continuous(name=orig.names[1])+
    scale_y_continuous(name=orig.names[2])+
    guides(alpha = FALSE, size = FALSE)+ theme_bw()
} 



# ## find the bounds for the complete data 
# xbnds <- range(c(A$x, B$x))
# ybnds <- range(c(A$y, B$y))
# nbins <- 30
# 
# #  function to make a data.frame for geom_hex that can be used with stat_identity
# 
# 
# library(ggplot2)
# #want to compute the histogram estimating the difference of the two densities
# ggplot(df, aes(x1,x2, group=group)) + geom_hex(bins=5, alpha=.2)# + facet_wrap(~group)
# 
# library(dplyr)
# p1 <- filter(df, group==1) %>%
#   ggplot(aes(x1,x2)) + geom_hex(bins=7)
# 
# 
# 
# library(hexbin)
# library(ggplot2)
# 
# set.seed(2)
# xA <- rnorm(1000,1)
# set.seed(3)
# yA <- rnorm(1000,1)
# set.seed(4)
# hbinA <- hexbin(xA, yA, xbins = 40, IDs = TRUE)
# 
# A <- data.frame(x = xA, y = yA)
# 
# set.seed(5)
# xB <- rnorm(1000)
# set.seed(6)
# yB <- rnorm(1000)
# set.seed(7)
# zB <- sample(c(1, 0), 20, replace = TRUE, prob = c(0.4, 0.6))
# hbinB <- hexbin(xB, yB, xbins = 40, IDs = TRUE)
# 
# B <- data.frame(x = xB, y = yB)
# 
# 
# ggplot(A, aes(x, y, z = z)) + stat_summary_hex(fun = function(z) sum(z)/length(z), alpha = 0.8) +
#   scale_fill_gradientn(colours = c("blue","red")) +
#   guides(alpha = FALSE, size = FALSE)
# 
# ggplot(B, aes(x, y, z = z)) + stat_summary_hex(fun = function(z) sum(z)/length(z), alpha = 0.8) +
#   scale_fill_gradientn (colours = c("blue","red")) +
#   guides(alpha = FALSE, size = FALSE)
# 
# ## find the bounds for the complete data 
# xbnds <- range(c(A$x, B$x))
# ybnds <- range(c(A$y, B$y))
# nbins <- 30
# 
# #  function to make a data.frame for geom_hex that can be used with stat_identity
# makeHexData <- function(df) {
#   h <- hexbin(df$x, df$y, nbins, xbnds = xbnds, ybnds = ybnds, IDs = TRUE)
#   data.frame(hcell2xy(h),
#              z = h@count,
#              cid = h@cell)
# }
# 
# Ahex <- makeHexData(A)
# Bhex <- makeHexData(B)
# 
# ##  not all cells are present in each binning, we need to merge by cellID
# byCell <- merge(Ahex, Bhex, by = "cid", all = T)
# 
# ##  when calculating the difference empty cells should count as 0
# byCell$z.x[is.na(byCell$z.x)] <- 0
# byCell$z.y[is.na(byCell$z.y)] <- 0
# 
# ##  make a "difference" data.frame
# Diff <- data.frame(x = ifelse(is.na(byCell$x.x), byCell$x.y, byCell$x.x),
#                    y = ifelse(is.na(byCell$y.x), byCell$y.y, byCell$y.x),
#                    z = byCell$z.x - byCell$z.y)
# 
# ##  plot the results
# 
# ggplot(Ahex) +
#   geom_hex(aes(x = x, y = y, fill = z),
#            stat = "identity", alpha = 0.8) +
#   scale_fill_gradientn (colours = c("blue","red")) +
#   guides(alpha = FALSE, size = FALSE)
# 
# ggplot(Bhex) +
#   geom_hex(aes(x = x, y = y, fill = z),
#            stat = "identity", alpha = 0.8) +
#   scale_fill_gradientn (colours = c("blue","red"), trans="log") +
#   guides(alpha = FALSE, size = FALSE)
# 
# library(scales)
# tanh_trans <- function(scale){
#   trans <- function(x) tanh
#   inv <- atanh
#   scales::trans_new("tanh", transform = trans, inverse = inv)
# }
# 
# ggplot(Diff) +
#   geom_hex(aes(x = x, y = y, fill = z),
#            stat = "identity", alpha = 0.8) +
#   scale_fill_gradient2(low="blue",high="red", mid="white",
#                         trans=scales::trans_new("tanh", 
#                                                 transform = function(x) tanh(x/10), 
#                                                 inverse = function(x) atanh(x)*10),
#                         breaks=c(-10,-5,0,5,10)) +
#   guides(alpha = FALSE, size = FALSE)+ theme_bw()
