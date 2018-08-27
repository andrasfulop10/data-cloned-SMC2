clear all;

Eh=.025^2;

alpha1=.7;
alpha2=.1;
theta=1;

alpha0=Eh.*(1-alpha1-alpha2*(1+theta^2));

mu=.1/252;

beta=1;

%simulate data
T=1000;

hT=zeros(T,1);
lnP_intrinsic=zeros(T,1);
P=zeros(T,1);

observationpattern=ones(T,1);
observationpattern(100:103)=0;

observationpattern(500:506)=0;

hcurrent=Eh;

hT(1)=hcurrent;
lnP_intrinsic(1)=log(100);
P(1)=exp(lnP_intrinsic(1));

B=.05;

for t=2:T
   
    epsilon=GEDSTD_invcdf(rand,beta);
    
    rt=mu+sqrt(hcurrent)*epsilon;
    
    hcurrent=alpha0+alpha1*hcurrent+alpha2*hcurrent*(epsilon-theta)^2;
    
    hT(t)=hcurrent;
    
    lnP_intrinsic(t)=lnP_intrinsic(t-1)+rt;
    
    if observationpattern(t)
    
        if observationpattern(t-1)
            
            LB=P(t-1)*(1-B);
            UB=P(t-1)*(1+B);
            
            P(t)=min(UB,max(LB,exp(lnP_intrinsic(t))));
            
        else
           
            P(t)=exp(lnP_intrinsic(t));
            
        end    
        
    else
    
        P(t)=nan;
        
    end
        
end

timestamp=(1:T)';

keep=find(isfinite(P));

hT=hT(keep);
P=P(keep);
lnP_intrinsic=lnP_intrinsic(keep);

timestamp=timestamp(keep);

Nparticles=50;

Nrepgap=20;

gpu_switch=1;

Nparam=1000;

Ti=50;

mu=repmat(mu,Nparam,1);
alpha0=repmat(alpha0,Nparam,1);
alpha1=repmat(alpha1,Nparam,1);
alpha2=repmat(alpha2,Nparam,1);
theta=repmat(theta,Nparam,1);
beta=repmat(beta,Nparam,1);

tic

[l,FiltMean_h,FiltMean_lnPintrinsic]=...
    PF_CensoredGARCH_GED(P,timestamp,alpha0,alpha1,alpha2,theta,mu,beta,...
    B,Nparticles,gpu_switch,Ti,Nrepgap);

toc

