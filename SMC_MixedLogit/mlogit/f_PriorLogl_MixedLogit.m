function lnpdf_prior=f_PriorLogl_MixedLogit(X,prior_param)
% evaluate initialization samplers likelihood
% use normal initialization alpha,delta and truncated normal for sigmaV
% 4/24/16 prior_param.sig.^2 for f_normlnpdf
% 5/13/16 non-negativity 4:6

Nparam=size(X,1);

lnpdf_prior=sum(f_normlnpdf(X,repmat(prior_param.mu,Nparam,1),repmat(prior_param.sig.^2,Nparam,1)),2);

drop=find(X(:,4)<0 | X(:,5)<0 | X(:,6)<0);
 
lnpdf_prior(drop)=-inf;  
