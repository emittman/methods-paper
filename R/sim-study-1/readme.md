## Plan for simulation study 1

### Construction of simulated data:

1. Select random draw $\mathcal{P}^{(s)}$ from posterior distribution of $\mathcal{P}$

2. Sample $(\beta_g,\sigma^2_g)$ from $\mathcal{P}^{(s)}$

3. Sample $y_{g,rep} \sim N(X\beta_g,\sigma^2_g)$

Produce 6-10 such data sets.

### Analysis

Run Gibbs sampler for each data set.

1. look at coverage of posterior credible intervals for gene-specific parameters

2. look at calibration of posterior probabilities for proposition(s) of interest

3. compare MSE to unpooled/unbiased estimateors

4. compare ROC curves to limma