function [l,filtersettings,FiltMeans]=Get_logl_GARCH(X,Data,filtersettings)

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

%logreturns
Y=log(Data.P(2:end)./Data.P(1:end-1));

[li,FiltMeans]=...
    Filter_GARCH(Y,alpha0,alpha1,alpha2,theta,mu,beta,filtersettings.gpu_switch);

l=sum(li(:,filtersettings.burnin+1:end),2);




