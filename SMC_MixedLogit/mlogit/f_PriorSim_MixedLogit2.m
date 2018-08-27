function initial_sample=f_PriorSim_MixedLogit2(prior_param,Nparam)
% alpha delta sigmaV
%
% use normal initialization alpha and truncated normal for delta and sigmaV
% 2/28/16 for mixed logit
% 5/13/16 lognormal 4:6

%A=transpose(sqrtm(prior_param.sig));
%initial_sample=repmat(prior_param.mu,Nparam,1)+randn(Nparam,4)*A;

initial_sample=mvnrnd(prior_param.mu,prior_param.sig,Nparam);
initial_sample(:,3:4)=exp(initial_sample(:,3:4));
