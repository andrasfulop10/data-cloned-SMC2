function [l,filtersettings]=KF_logl(X,Y,filtersettings)
%[mu phi log(sigma_eta) log(sigma_epsilon)]

mu=X(:,1);

phi=X(:,2);

%variance of state noise
sigma2_eta=exp(2*X(:,3));

%variance of measurmemnt noise
sigma2_epsilon=exp(2*X(:,4));

%initialize hidden state at stationary distribution

E_alpha=mu./(1-phi);
sig2_alpha=sigma2_eta./(1-phi.^2);


[li,Particles.E_alpha,Particles.sig2_alpha]=...
    KF(mu,phi,sigma2_eta,sigma2_epsilon,Y,E_alpha,sig2_alpha);


l=sum(li,2);




