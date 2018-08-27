function [ESS,normfac,AcceptRate,X,l,lnw,lnpdf_prior,gamma,filtersettings,smcsettings]=...
    SMC_Step(X,l,lnw,lnpdf_prior,prior_param,gamma,...
    Y,filtersettings,smcsettings)

%one step of the density-tempered SMC
AcceptRate=nan;

normfac=0;

    gammavec=linspace(-100,0,1000);
    
    ESSvec=TempFun(gammavec,lnw,l-lnpdf_prior);
    
    j=find(ESSvec<smcsettings.ESS_bound);
    
    if ~isempty(j)
        
        gammaincrement=gammavec(j(1));
        
        gammanew=min(gamma+exp(gammaincrement),1);
        
    elseif min(ESSvec)>smcsettings.ESS_bound
        
        gammanew=1;
        
    end
            
    incrementalw=(gammanew-gamma)*(l-lnpdf_prior);
    
    j=find(isnan(incrementalw) | imag(incrementalw)~=0);
    
    incrementalw(j)=-inf;
    
    %get incremental normalizing ratio
    W_prev=exp(lnw-max(lnw)); W_prev=W_prev/sum(W_prev);
    
    max_incrementalw=max(incrementalw);
    
    normfac=normfac+log(sum(W_prev.*exp(incrementalw-max_incrementalw)))+max_incrementalw;
    %end of normalizing ratio computations
    
    lnw=lnw+incrementalw;
    
    %compute normalized weights and compute ESS
    W=exp(lnw-max(lnw));
    
    W=W/sum(W);
    
    ESS=1/sum(W.^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %resample and move if ESS<ESS_bound%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ESS<smcsettings.ESS_bound
        
        if smcsettings.verbose
        
            disp(['Resample at gamma=' num2str(gammanew)]);
        end
        %%%%%%%%%%
        %Resample%
        %%%%%%%%%%
        
        [X,lnw,l,lnpdf_prior]=ResampleSet(X,lnw,l,lnpdf_prior);
               
        counter=0;
        
        while (counter<smcsettings.Nmoves)
            
            counter=counter+1;
            
            [X,l,lnpdf_prior,AcceptRate,filtersettings,smcsettings]=...
                MoveSet(X,prior_param,lnpdf_prior,l,...
                Y,gammanew,filtersettings,smcsettings);
                   
            if smcsettings.verbose
            
            %disp(['unique x after move # ' num2str(counter)]);
            
            %size(unique(X,'rows'),1)
                      
            end
        end
            
        
    end
    
    
    gamma=gammanew;
  
    





