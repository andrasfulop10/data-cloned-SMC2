function prior_param=PriorFit(X)

prior_param.mu=mean(X);

prior_param.V=cov(X);
