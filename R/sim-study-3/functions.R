run_mcmc <- function(data, design){
  require(cudarpackage)
  require(limma)
  require(edgeR)
  require(dplyr)
  # dge <- DGEList(data) %>% calcNormFactors(method="TMM")
  # eff.lib.size <- dge[[2]]$lib.size * dge[[2]]$norm.factors
  G <- nrow(data)
  log_cpm <- log(data +1) #- log2(matrix(rep(eff.lib.size, each=G), G, 16)+1) + log2(1e6)

  dat <- formatData(counts=log_cpm, X=design, transform_y=identity, voom = FALSE, normalize = FALSE)
  ind_est <- indEstimates(dat)
  priors <- formatPriors(K=2^10, estimates=ind_est, A=3, B=3/sqrt(dat$G))
  C <- list(b2 = matrix(c(0,1,0,0,0), 1,5,byrow=T),
            b3 = matrix(c(0,0,1,0,0), 1,5,byrow=T),
            b4 = matrix(c(0,0,0,1,0), 1,5,byrow=T),
            b5 = matrix(c(0,0,0,0,1), 1,5,byrow=T)
  )
  contr <- formatControl(n_iter = 2,
                         thin = 1,
                         warmup = 30000,
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
  
  contr$n_iter <- as.integer(100000)
  contr$thin <- as.integer(50)
  contr$warmup <- as.integer(50000)
  contr$idx_save <- 0:(dat$G-1)
  contr$n_save_P <- as.integer(100)
  contr$methodPi <- "stickBreaking"
  
  samples <- mcmc(dat, priors, contr, init_chain)
  samples
}
