initialize_chain <- function(seed, methodPi, n.init, n.warmup, n.sample){

  cat(methodPi)  
    
  set.seed(seed)
  load("data/cuda-dat-heterosis.Rdata")
  load("data/ind-est-heterosis.RData")
  priors <- formatPriors(K=2^12*1.5, estimates = ind_est, A=5, B=5/sqrt(cuda_dat$G))
  
  C <- list(extreme_heterosis = matrix(c(0, 1, 1, 1, 0,
			 	 0, 1, 1,-1, 0,
				 0,-1, 1, 1, 0,
				 0,-1, 1,-1, 0),4, 5, byrow=T)) 
  
  contr <- formatControl(n_iter = 2,
                         thin = 1,
                         warmup = n.init,
                         methodPi = "symmDirichlet",
                         idx_save = 1,
                         n_save_P = 1,
                         alpha_fixed = FALSE)
  
  #run a pilot chain and reorder clusters
  start.chain <- initFixedGrid(priors = priors, estimates = ind_est, C = C)
  init.run <- mcmc(cuda_dat, priors, contr, start.chain)
  
  id <- order(init.run[['state']]$pi, decreasing=TRUE)
  init.chain <- with(init.run[['state']],
                     formatChain(beta = beta[,id],
                                 pi = exp(pi[id]),
                                 tau2 = tau2[id],
                                 zeta =  start.chain$zeta,
                                 alpha = alpha,
                                 C = C))
  
  contr$n_iter <- as.integer(n.sample)
  contr$thin <- as.integer(100)
  contr$warmup <- as.integer(n.warmup)
  contr$idx_save <- as.integer(0:(cuda_dat$G-1))
  contr$n_save_P <- as.integer(100)
  contr$methodPi <- methodPi
  
  list(priors = priors,
       control = contr,
       seed = .Random.seed,
       init.chain = init.chain)
}

  
sample_bnp_model <- function(settings){
  load("data/cuda-dat-heterosis.Rdata")
  set.seed(settings$seed)
  samples <- with(settings, mcmc(cuda_dat, priors, control, init.chain))
  samples
}
