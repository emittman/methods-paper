#' @name sim_study_2_setup
#' @return list containing estimates from voom-limma pipeline to use as truth in simulation

sim_study_2_setup <- function(){
  load("../data/heterosis_counts.RData")
  load("../data/heterosis_design.RData")
  require(limma)
  require(edgeR)
  require(dplyr)
  
  dge <- DGEList(counts = my_dat_wide[,2:17]) %>%
    calcNormFactors(method = "TMM")
  
  voom_dge <- voom(counts = dge, design = X, normalize.method = "none")
  fit <- lmFit(voom_dge)
  fit_eb <- eBayes(fit)
  
  prereqs <- list(betas = fit$coefficients,
                  sigma2 = fit_eb$s2.post,
                  weights = voom_dge$weights,
                  eff.lib.size = with(dge$samples, lib.size*norm.factors))
  prereqs
}

#' @title generate_log_counts
#' @par oldG number of genes in model fit provided in prereqs
#' @par newG number of genes in simulated data set
#' @par design N by p matrix
#' @par prereqs list containing true values, from sim_study2_setup
#' @return list containing simulated counts and a list containing the true parameter values 

generate_counts <- function(oldG, newG, design, prereqs){
  N <- nrow(design)
  samp_id <- sample(oldG, newG*1.5)
  
  beta_rep <- prereqs$beta[samp_id,]
  sigma2_rep <- prereqs$sigma2[samp_id]
  weights_rep <- prereqs$weights[samp_id,]
  offset <- mean(log2(prereqs$eff.lib.size)) - log2(1e6)
  ystar <- t(sapply(1:(newG*1.5), function(g){
    rnorm(16, X%*%beta_rep[g,] + offset, sqrt(sigma2_rep[g]/weights_rep[g,]))
  }))
  y <- round(2^(ystar))
  zero_id <- which(apply(y, 1, function(g) all(g==0)))
  stopifnot(length(zero_id) <= newG*1.25)
  if(length(zero_id)>0){
    y <- y[-zero_id,]
    beta_rep <- beta_rep[-zero_id]
    sigma2_rep <- sigma2_rep[-zero_id]
    weights_rep <- weights_rep[-zero_id,]
  }
  beta <- beta_rep[1:newG,]
  sigma2 <- sigma2_rep[1:newG]
  weights <- weights_rep[1:newG,]
  y <- y[1:newG,]
  sim <- list(y=y,
              truth = list(beta = beta,
                           sigma2 = sigma2,
                           weights = weights))
  sim
}

#'@title estGammaPrior
#'@par x vector of sample variances
#'@return named vector: a = shape and b = scale
estGammaPrior <- function(x){
  a <- 1
  b <- 1
  for(i in 1:1000){
    anew <- mean(1/x+.0001, tr=.1)*b
    bnew <- 1/( quantile(1/x, probs=.95)/qgamma(.95,anew,1) )
    if(abs(anew-a)<.0001 & abs(bnew-b)<.0001){
      break
    }
    a <- anew
    b <- bnew
  }
  return(c(a=a,b=b))
}

#' @title run_mcmc
#' @par data
#' @return myMcmcObj

run_mcmc <- function(data, design, voom=TRUE){
  require(cudarpackage)
  require(limma)
  require(edgeR)
  require(DESeq2)
  #use deseq to estimate priors
  deseq.dat <- DESeqDataSetFromMatrix(data$y, colData = as.data.frame(design[,-1]),
                                      design = ~ parent_hd+hybrid+hybrid_hd+flow_cell)
  fit.deseq <- DESeq(deseq.dat)
  mles <- estimateMLEForBetaPriorVar(fit.deseq)
  priors.beta <- estimateBetaPriorVar(mles, betaPriorMethod = "quantile")
  
  if(!voom){
    dge <- DGEList(data$y) %>% calcNormFactors(method="TMM")
    eff.lib.size <- dge[[2]]$lib.size * dge[[2]]$norm.factors
    G <- nrow(data$y)
    log_cpm <- log2(data$y+.5) - log2(matrix(rep(eff.lib.size, each=G), G, 16)+1) + log2(1e6)
    data$y <- log_cpm
  }
  
  dat <- formatData(counts=data$y, X=design, transform_y=identity, voom = voom, normalize = voom)
  ind_est <- indEstimates(dat)
  priors.sigma2 <- estGammaPrior(ind_est$sigma2)
  
  priors <- formatPriors(K=2^12, estimates=ind_est, A=3, B=3/sqrt(dat$G))
  priors$lambda2[2:5] <- 1/priors.beta[2:5]
  priors$a <- priors.sigma2[1]
  priors$b <- priors.sigma2[2]
  
  C <- list(high_mean = matrix(c(0,1,1,0,0,
                                 0,-1,1,0,0), 2,5,byrow=T),
            hp_h12    = matrix(c(0,1,1,1,0,
                                 0,-1,1,1,0), 2,5,byrow=T),
            lp_h12    = matrix(c(0,1,-1,-1,0,
                                 0,-1,-1,-1,0), 2, 5, byrow=T),
            hp_h21    = matrix(c(0,1,1,-1,0,
                                 0,-1,1,-1,0), 2,5, byrow=T),
            lp_h21    = matrix(c(0,1,-1,1,0,
                                 0,-1,-1,1,0), 2, 5, byrow=T),
            de_p1     = matrix(c(0,1,0,0,0), 1, 5, byrow=T),
            high_or_mid = matrix(c(0,0,1,0,0), 1, 5, byrow=T)
  )
  contr <- formatControl(n_iter = 2,
                         thin = 1,
                         warmup = 10000,
                         methodPi = "symmDirichlet",
                         idx_save = 1,
                         n_save_P = 1,
                         alpha_fixed = FALSE)
  
  start_chain <- initFixedGrid(priors = priors, estimates = ind_est, C = C)
  init_run <- mcmc(dat, priors, contr, start_chain)
  
  id <- order(init_run[['state']]$pi, decreasing=TRUE)
  init_chain <- formatChain(beta = init_run[['state']]$beta[,id],
                            pi = exp(init_run[['state']]$pi[id]),
                            tau2 = init_run[['state']]$tau2[id],
                            zeta =  start_chain$zeta,
                            alpha = init_run[['state']]$alpha,
                            C = C)
  
  contr$n_iter <- as.integer(50000)
  contr$thin <- as.integer(25)
  contr$warmup <- as.integer(10000)
  contr$idx_save <- 0:(dat$G-1)
  contr$n_save_P <- as.integer(100)
  contr$methodPi <- "stickBreaking"
  
  samples <- mcmc(dat, priors, contr, init_chain)
  samples
}

