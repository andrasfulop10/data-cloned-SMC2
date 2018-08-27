function X=PriorSim(prior_param,Nparam)
% mu phi log(sigma_eta) log(sigma_epsilon)

X=mvnrnd(prior_param.mu,prior_param.V,Nparam*10);

keep=find(abs(X(:,2))<1);

X=X(keep(1:Nparam),:);
