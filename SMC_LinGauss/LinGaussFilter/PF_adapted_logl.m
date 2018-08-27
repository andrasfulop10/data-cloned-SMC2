function [l,filtersettings]=PF_adapted_logl(X,Y,filtersettings)

if filtersettings.gpu_switch
    X=gpuArray(X);
end

mu=X(:,1);

Nparam=length(mu);

phi=X(:,2);

%variance of state noise
sigma2_eta=exp(2*X(:,3));

%variance of measurmemnt noise
sigma2_epsilon=exp(2*X(:,4));

%initialize hidden state at stationary distribution

E_alpha=mu./(1-phi);
sig_alpha=sqrt(sigma2_eta./(1-phi.^2));

if filtersettings.randomnumbers.CRN==1
    
    randnorm=filtersettings.randomnumbers.epsilon(:,1);
    
    if filtersettings.gpu_switch
        
        randnorm=gpuArray(randnorm);
    end
    
    randnorm=repmat(randnorm,1,Nparam);
    
else
    
    if filtersettings.gpu_switch
        randnorm=gpuArray.randn(filtersettings.Nparticles,Nparam);
    else
        randnorm=randn(filtersettings.Nparticles,Nparam);
    end
    
end

alpha0=repmat(E_alpha',filtersettings.Nparticles,1)+repmat(sig_alpha',filtersettings.Nparticles,1).*randnorm;

li=PF_adapted(mu,phi,sigma2_eta,sigma2_epsilon,Y,alpha0,filtersettings);

l=sum(li,2);




