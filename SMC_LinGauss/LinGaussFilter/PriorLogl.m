function lnpdf_prior=PriorLogl(X,prior_param)
% mu phi log(sigma_eta) log(sigma_epsilon)

lnpdf_prior=log(mvnpdf(X,prior_param.mu,prior_param.V));

drop=find(X(:,2)<-1 | X(:,2)>1);

lnpdf_prior(drop)=-inf;
