#' @title generate_log_counts
#' @par G number of genes
#' @par design N by p matrix
#' @par P.samples set of draws from the posterior for P (myMcmcObj$samples$P)
#' @return list containing simulated log-counts and a list containing the true parameter values 

generate_log_counts <- function(G, design, P.samples){
  N <- nrow(design)
  
  # if(!is.null(weights)) stopifnot(G*N == length(weights))
  
  P.dims <- dim(P.samples)
  K <- P.dims[1]
  p <- P.dims[2]-2
  S <- P.dims[3]
  
  #random draw of P
  id <- sample(S, 1)
  P.truth <- P.samples[,,id]
  
  #sample zeta
  zeta <- sample(K, G, replace=TRUE, prob = P.truth[,1])
  beta.truth <- P.truth[zeta, 2:(p+1)]
  sigma.truth <- 1/sqrt(P.truth[zeta, p+2])
  # if(is.null(weights)){
    y <- t(sapply(1:G, function(g) rnorm(N, design %*% beta.truth[g,], sigma.truth[g])))
  # }
  list(y=y, truth=list(beta = beta.truth, sigma = sigma.truth))
}


#' @title run_mcmc
#' @par data
#' @return myMcmcObj

run_mcmc <- function(data, design){
  require(cudarpackage)
  dat <- formatData(counts=data$y, X=design, transform_y=identity, voom = FALSE)
  ind_est <- indEstimates(dat)
  priors <- formatPriors(K=2^12, estimates=ind_est, A=3, B=3/sqrt(dat$G))
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
            de_p1     = matrix(c(0,1,0,0,0), 1, 5, byrow=T)
            )
  contr <- formatControl(n_iter = 2,
                         thin = 1,
                         warmup = 5000,
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
  
  contr$n_iter <- as.integer(30000)
  contr$thin <- as.integer(15)
  contr$warmup <- as.integer(10000)
  contr$idx_save <- 0:(dat$G-1)
  contr$n_save_P <- as.integer(100)
  contr$methodPi <- "stickBreaking"
  
  samples <- mcmc(dat, priors, contr, init_chain)
  samples
}

#fit_limma <- 

