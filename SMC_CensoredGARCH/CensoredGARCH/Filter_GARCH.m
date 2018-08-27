function [l,hstore]=...
    Filter_GARCH(Y,alpha0,alpha1,alpha2,theta,mu,beta,gpu_switch)
%function runs particle filter for censored NGARCH(1,1) model
%Y are logreturns
if gpu_switch
    alpha0=gpuArray(alpha0);
    alpha1=gpuArray(alpha1);
    alpha2=gpuArray(alpha2);
    theta=gpuArray(theta);
    mu=gpuArray(mu);
    beta=gpuArray(beta);
end

T=length(Y);

Nparam=length(alpha0);

l=zeros(Nparam,T);

hstore=zeros(Nparam,T);

%initialize the states from variance

%h=var(log(1+Ret)).*(0.5+rand(Nparticles,Nparam));

%initialize states from unconditional mean
h=alpha0./(1-alpha1-alpha2.*(1+theta.^2));

%run time loop with bootstrap filter
for t=1:T
         
    epsilon=(Y(t)-mu)./sqrt(h);
    
    l(:,t)=GEDSTD_lnpdf(epsilon,beta)-log(h)/2;
               
    %GARCH equation
    h=alpha0+alpha1.*h+alpha2.*h.*(epsilon-theta).^2;
        
    hstore(:,t)=h;
    
end

if gpu_switch
    
    l=gather(l);
    hstore=gather(hstore);
end





