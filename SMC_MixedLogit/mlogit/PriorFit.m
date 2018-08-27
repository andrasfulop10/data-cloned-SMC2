function prior_param=PriorFit(X)
X(:,3:4)=log(X(:,3:4));
prior_param.mu=mean(X);

prior_param.sig=cov(X);
