function l=PF_adapted(mu,phi,sigma2_eta,sigma2_epsilon,Y,alpha0,filtersettings)

[Nparticles,Nparam]=size(alpha0);

T=length(Y);

mu=repmat(mu',Nparticles,1);

phi=repmat(phi',Nparticles,1);

sigma2_eta=repmat(sigma2_eta',Nparticles,1);

sigma2_epsilon=repmat(sigma2_epsilon',Nparticles,1);

l=zeros(Nparam,T);

if filtersettings.gpu_switch
    lnw_equal=gpuArray.ones(Nparticles,Nparam)*log(1/Nparticles);
else    
    lnw_equal=ones(Nparticles,Nparam)*log(1/Nparticles);
end

lnw=lnw_equal;

alpha_current=alpha0;

for t=1:T

    
    %%%%%%%%%%%%%%%%%%%%%
    %forecasted moments
    %%%%%%%%%%%%%%%%%%%%%
    E_newalpha=mu+phi.*alpha_current;
    
    V_newalpha=sigma2_eta;
    
    V_Y=sigma2_eta+sigma2_epsilon;
    
    Cov_stateY=sigma2_eta;
    
    E_Y=E_newalpha;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %include info from new observation%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lnw=lnw+normlnpdf(Y(t),E_Y,V_Y);
      
    %compute likelihood and resample
    maxlnw=max(lnw);

    lnw=bsxfun(@minus,lnw,maxlnw);

    w=exp(lnw);

    w(not(isfinite(w)))=0;

    cumweights=cumsum(w);

    sumw=cumweights(end,:);

    logl_t=(log(sumw)+maxlnw)';
    
    if filtersettings.gpu_switch
        logl_t=gather(logl_t);
    end
    
    l(:,t)=logl_t;
    
    cumweights=bsxfun(@rdivide,cumweights,sumw);
    
    if filtersettings.randomnumbers.CRN==1
        
        randu=filtersettings.randomnumbers.u(:,t);
        
        if filtersettings.gpu_switch
            
            randu=gpuArray(randu);
        end
        
        randu=repmat(randu,1,Nparam);
        
    else
        
        if filtersettings.gpu_switch
            randu=gpuArray.rand(Nparticles,Nparam);
        else
            randu=rand(Nparticles,Nparam);
        end
    end
        
    v=(0:Nparticles-1)';
    
    if filtersettings.gpu_switch
        v=gpuArray(v);
    end
    
    u=(repmat(v,1,Nparam)+randu)/Nparticles;
    
    bin=Resample_Mtx(u,cumweights);
    
    mtx_bin=repmat((0:Nparam-1)*Nparticles,Nparticles,1)+bin;
    
    E_newalpha=E_newalpha(mtx_bin);
    
    V_newalpha=V_newalpha(mtx_bin);
    
    V_Y=V_Y(mtx_bin);
    
    Cov_stateY=Cov_stateY(mtx_bin);
    
    E_Y=E_Y(mtx_bin);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now update conditional state moments and sample from new states 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_newalpha=E_newalpha+Cov_stateY./V_Y.*(Y(t)-E_Y);
    
    V_newalpha=V_newalpha-Cov_stateY.^2./V_Y;
    
     if filtersettings.randomnumbers.CRN==1
        
        randnorm=filtersettings.randomnumbers.epsilon(:,t+1);
        
        if filtersettings.gpu_switch
            
            randnorm=gpuArray(randnorm);
        end
        
        randnorm=repmat(randnorm,1,Nparam);
        
    else
                
        if filtersettings.gpu_switch
            randnorm=gpuArray.randn(Nparticles,Nparam);
        else
            randnorm=randn(Nparticles,Nparam);
        end
    end
        
    alpha_current=E_newalpha+sqrt(V_newalpha).*randnorm;
    
    lnw=lnw_equal;
    
end



