function [l,FiltMean_h,FiltMean_lnPintrinsic]=...
    PF_CensoredGARCH_GED(P,timestamp,alpha0,alpha1,alpha2,theta,mu,beta,...
    Bounds,Nparticles,gpu_switch,Ti,Nrepgap)
%function runs particle filter for censored NGARCH(1,1) model

if gpu_switch
    alpha0=gpuArray(alpha0);
    alpha1=gpuArray(alpha1);
    alpha2=gpuArray(alpha2);
    theta=gpuArray(theta);
    mu=gpuArray(mu);
    beta=gpuArray(beta);
    
end

T=length(P);

Nparam=length(alpha0);

alpha0=repmat(alpha0',Nparticles,1);
alpha1=repmat(alpha1',Nparticles,1);
alpha2=repmat(alpha2',Nparticles,1);
theta=repmat(theta',Nparticles,1);
mu=repmat(mu',Nparticles,1);
beta=repmat(beta',Nparticles,1);

l=zeros(Nparam,T-1);

FiltMean_h=zeros(Nparam,T-1);

FiltMean_lnPintrinsic=zeros(Nparam,T-1);

%initialize the states from variance

%initialize states from unconditional mean and then iterate for Ti periods
%to arrive to stationary distribution of variance
h=alpha0./(1-alpha1-alpha2.*(1+theta.^2));

for t0=1:Ti
    
    if gpu_switch
        u=gpuArray.rand(Nparticles,Nparam);
    else
        u=rand(Nparticles,Nparam);
    end
        
    epsilon=GEDSTD_invcdf(u,beta);
    
    h=alpha0+alpha1.*h+alpha2.*h.*(epsilon-theta).^2;
end

%initialize intrinsic price process at first observed value
if gpu_switch
     lnP_intrinsic=log(P(1))*gpuArray.ones(Nparticles,Nparam);
else
     lnP_intrinsic=log(P(1))*ones(Nparticles,Nparam);
end

%run time loop 
for t=2:T
    
    timegap=timestamp(t)-timestamp(t-1);
    
    if timegap==1
        loglik_y=PredictiveLikl_step1(P(t),P(t-1),mu,beta,Bounds,h,lnP_intrinsic);
    
        lnw=loglik_y-log(Nparticles);
    
        [l(:,t),mtx_bin]=Resample(lnw,gpu_switch,Nparticles);
        
        h=h(mtx_bin);
                  
        lnP_intrinsic=lnP_intrinsic(mtx_bin);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %propagate states given the past and the current observation%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [h,lnP_intrinsic]=PropagateVol(h,lnP_intrinsic,P(t),P(t-1),alpha0,alpha1,alpha2,theta,mu,beta,...
                Bounds,gpu_switch);
                        
    else
        %simulate forward Nrepgap*Nparticles h, P up to before the new time
        %points and compute new likelihood
        %loglik_y,P_beforet,h_beforet are all (Nparticles*Nrepgap)*Nparam
        %matrices
        [loglik_y,lnP_intrinsic_beforet,h_beforet]=...
            PredictiveLikl_stepk(P(t),alpha0,alpha1,alpha2,theta,mu,beta,...
                Bounds,h,lnP_intrinsic,Nrepgap,timegap,gpu_switch);
                    
        lnw=loglik_y-log(Nparticles*Nrepgap);
           
        
        [l(:,t),mtx_bin]=Resample(lnw,gpu_switch,Nparticles);
                
        h_beforet=h_beforet(mtx_bin);   
        lnP_intrinsic_beforet=lnP_intrinsic_beforet(mtx_bin);
                    
        [h,lnP_intrinsic]=PropagateVol_stepk(P(t),lnP_intrinsic_beforet,h_beforet,...
        alpha0,alpha1,alpha2,theta,mu,gpu_switch) ;   
             
    end
    
       mean_h=mean(h);
       
       mean_lnPintrinsic=mean(lnP_intrinsic);
    
       if gpu_switch
         mean_h=gather(mean_h);
         
         mean_lnPintrinsic=gather(mean_lnPintrinsic);
       end
    
    FiltMean_h(:,t-1)=mean_h';
    FiltMean_lnPintrinsic(:,t-1)=mean_lnPintrinsic';
        
end





