#Simulation study 2

- Use posterior means as truth, reverse engineer log_counts conditioned on voom weights,
i.e. $y_{gi,rep} \sim N\left( \operatorname{E}(\mu_{gi}|y), \operatorname{E}(\sigma_{g}|y)^2/w_{gi})\right)$

- Exponentiate (base 2) and round to get counts

- Compare ROC and MSE for various methods

- Issues to resolve: how to rank genes for other methods with composite null structure