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
  stopifnot(length(zero_id) <= newG*1.5)
  y <- y[-zero_id,]
  y <- y[1:newG,]
  sim <- list(y=y,
              truth = list(beta = beta_rep,
                           sigma2 = sigma2_rep,
                           weights = weights_rep))
  sim
}


#' @title run_mcmc
#' @par data
#' @return myMcmcObj

run_mcmc <- function(data, design, voom=TRUE){
  require(cudarpackage)
  require(limma)
  require(edgeR)
  if(!voom){
    dge <- DGEList(data$y) %>% calcNormFactors(method="TMM")
    eff.lib.size <- dge[[2]]$lib.size * dge[[2]]$norm.factors
    G <- nrow(data$y)
    log_cpm <- log2(data$y+.5) - log2(matrix(rep(eff.lib.size, each=G), G, 16)+1) + log2(1e6) 
  }
  
  dat <- formatData(counts=data$y, X=design, transform_y=identity, voom = voom, normalize = voom)
  ind_est <- indEstimates(dat)
  priors <- formatPriors(K=2^11*1.5, estimates=ind_est, A=3, B=3/sqrt(dat$G))
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

