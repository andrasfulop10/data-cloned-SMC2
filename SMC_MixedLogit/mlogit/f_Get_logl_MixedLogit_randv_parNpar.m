function [logl,filtersettings]=...
    f_Get_logl_MixedLogit_randv_parNpar(X,dataX,dataY,filtersettings)
%likelihood function of mixed logit also parallelized in fixed param
%direction

[J,N,K]=size(dataX);

M=filtersettings.Nparticles;  % # of simulation
Npar=size(X,1);  % # of samples

gpu_switch=filtersettings.gpu_switch;

if gpu_switch
    %bring input data on gpu dataY is not needed as it is only used to
    %index into matrices
    X=gpuArray(X);
    dataX=gpuArray(dataX);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulate betas (N*K*M*Npar)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=permute(repmat(X(:,1:K),[1 1 M N]),[4 2 3 1]); %Npar*K to (N*K*M*Npar)
sigma=permute(repmat(X(:,K+1:end),[1 1 M N]),[4 2 3 1]); %Npar*K to (N*K*M*Npar)

if gpu_switch
    epsilon=gpuArray.randn(N,K,M,Npar);
else
    epsilon=randn(N,K,M,Npar);
end

beta=mu+sigma.*epsilon;

%the coefficient of first coefficient is lognormal
beta(:,1,:,:)=exp(beta(:,1,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute exponential of the cross-product of beta and dataX%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(N*K*M*Npar*J)
dataX=permute(repmat(dataX,[1 1 1 M Npar]),[2 3 4 5 1]);

%(N*M*Npar*J)
Xbeta=squeeze(sum(bsxfun(@times,dataX,beta),2));

%normalize Xbeta by the mean in J-direction. normalization has no
%effect on likelihood as we will divide with sum of exponentials in J dimension!
Xbeta=bsxfun(@minus,Xbeta,max(Xbeta,[],4));
        
expXbeta=exp(Xbeta);

%%%%%%%%%%%%%%%%%%%%%%%
%compute loglikelihood%
%%%%%%%%%%%%%%%%%%%%%%%

%compute normalizing factor N*M*Npar
sum_expXbeta=sum(expXbeta,4);

%unnormalized prob of alternative chosen
if gpu_switch
    expXbeta_I=gpuArray.zeros(N,M,Npar);
else
    expXbeta_I=zeros(N,M,Npar);
end

for i=1:N
   
    expXbeta_I(i,:,:)=expXbeta(i,:,:,dataY(i));
        
end

%N*Npar
lnPri=squeeze(log(mean(expXbeta_I./sum_expXbeta,2)));
        
%Npar*1
logl=sum(lnPri,1)';

if gpu_switch
    logl=gather(logl);
end






