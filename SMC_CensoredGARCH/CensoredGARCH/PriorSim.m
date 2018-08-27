function X=PriorSim(prior_param,Nparam)

X=mvnrnd(prior_param.mu,prior_param.V,Nparam);



