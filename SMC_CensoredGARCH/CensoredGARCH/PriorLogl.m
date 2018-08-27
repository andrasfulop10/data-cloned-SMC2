function lnpdf_prior=PriorLogl(X,prior_param)
%evaluate initialization samplers likelihood


Nparam=size(X,1);

lnpdf_prior=log(mvnpdf(X,prior_param.mu,prior_param.V));

