function lnpdf_prior=f_PriorLogl_MixedLogit2(X,prior_param)
% evaluate initialization samplers likelihood
% use normal initialization alpha,delta and truncated normal for sigmaV
% 4/24/16        prior_param.sig.^2 for f_normlnpdf
% 5/13/16 lognormal 4:6

Nparam=size(X,1);
    
X(:,3:4)=log(X(:,3:4));

lnpdf_prior=mvnlnpdf2(X,prior_param.mu,prior_param.sig);

logJacobian_prior=sum(-X(:,3:4),2);

lnpdf_prior=lnpdf_prior+logJacobian_prior;
    
    
    
    