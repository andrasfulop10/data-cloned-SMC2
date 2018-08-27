function logl=f_Get_logl_MixedLogit_fixv2(X,dataX,dataY,epsilon,Nparticles)
% Compute the likelihood of mixed logit
% fixed epsilon


[J,N,K]=size(dataX);  %J: alt, N:sample size, K:number of regressor

M=Nparticles;  %M: number of simulation
Npar=size(X,1);  % # of samples

logl=zeros(Npar,1);

dataX=permute(dataX,[2,3,1]);

%get indices of outcomes within expXbeta
C=(1:N)';

I1=repmat(C,1,M);

I2=repmat(dataY,1,M);

C=(1:M);

I3=repmat(C,N,1);

ind=sub2ind([N J M],I1,I2,I3);


for j=1:Npar

    mu=X(j,1:2);
    sigma=X(j,3:end); 
        

    %first M rows are betas for first individual etc...
    beta=repmat(mu,N*M,1)+epsilon.*repmat(sigma,N*M,1);
    beta(:,1)=exp(beta(:,1));%the first column correspond to lognormal

    beta=reshape(beta,[M N K]);

    %make beta N*K*M
    beta=permute(beta,[2 3 1]);

    %make beta N*K*J*M
    beta=permute(repmat(beta,[1 1 1 J]),[1 2 4 3]);

    %now compute cross-product of beta and X

    %N*J*M
    
    Xbeta=squeeze(sum(bsxfun(@times,dataX,beta),2));
    
    %normalize Xbeta by the mean in J-direction. normalization has no
    %effect on likelihood as we will divide the exponentials!
    Xbeta=bsxfun(@minus,Xbeta,max(Xbeta,[],2));
        
    expXbeta=exp(Xbeta);
    
    %N*M
    sum_expXbeta=squeeze(sum(expXbeta,2));

    %N*M
    expXbeta_I=expXbeta(ind);

    %N*1
    lnPr=log(mean(expXbeta_I./sum_expXbeta,2));
        
    logl(j)=sum(lnPr);
    
    
end

