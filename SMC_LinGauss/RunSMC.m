function [Xmean,Xsample,runtime,ESS,AcceptRate,Nparticles]=...
    RunSMC(Nparam,filtersettings,smcsettings,Y)
%Purpose: Execute one run of the cloned tempered SMC^2 routine

%simulate from initialization sampler and evaluate initialization density
X=feval(smcsettings.InitializationSim,smcsettings.Initialization_param,Nparam);

%allocate matrices to store results
lnw = zeros(Nparam,1);

l=zeros(Nparam,1);

Lparam=size(X,2);

cont=1;

i_k=1;

while cont
    
    tic
    %get to next clone number
    [ESS(i_k),AcceptRate(i_k),X,l,lnw,filtersettings,smcsettings]=...
        GetToNextK(X,l,i_k-1,lnw,Y,filtersettings,smcsettings);
    
    runtime(i_k)=toc;
    
    Nparticles(i_k)=filtersettings.Nparticles;
    
    W = exp(lnw - max(lnw));
    
    W = W/sum(W);
    
    Xmean(i_k,:)=sum(X.*repmat(W,1,Lparam));
    
    Xsample(:,:,i_k) = Resample_vec(X, W, Nparam);
        

    if (i_k==smcsettings.Kend)
             
        cont=0;
    
    end
        
    i_k=i_k+1;
    
end