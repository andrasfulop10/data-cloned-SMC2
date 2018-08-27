function [l,filtersettings,FiltMeans]=Get_logl(X,Data,filtersettings)

X=TransformParam(X);

if filtersettings.gpu_switch
    X=gpuArray(X);
end

alpha0=X(:,1);
alpha1=X(:,2);
alpha2=X(:,3);
theta=X(:,4);
mu=X(:,5);
beta=X(:,6);

[li,FiltMeans.h,FiltMeans.lnPintrinsic]=...
    PF_CensoredGARCH_GED(Data.P,Data.timestamp,alpha0,alpha1,alpha2,theta,mu,beta,...
    filtersettings.Bounds,...
    filtersettings.Nparticles,filtersettings.gpu_switch,filtersettings.Ti,...
    filtersettings.Nrepgap);

l=sum(li(:,filtersettings.burnin+1:end),2);




