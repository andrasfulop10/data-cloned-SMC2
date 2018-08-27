function [ESS,AcceptRate,X,l,lnw,filtersettings,smcsettings]=...
    f_GetToNextK(X,l,k_current,lnw,dataX,dataY,filtersettings,smcsettings)
%This function adds one more data cloning steps to sample

%number of resample-move steps
count_RM=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute loglikelihood estimate necessary for tempering up to next integer%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%obtain new estimate of loglikelihood
[l_nextk,filtersettings]=feval(filtersettings.filterfun,X,dataX,dataY,filtersettings); 

%if k_current=0, source density is a parametric initialization distribution
%with logl smcsettings.InitializationLogl 
if k_current==0
    
    initdensity_flag=1;
        
    lnpdf_init=feval(smcsettings.InitializationLogl,X,smcsettings.Initialization_param);
       
else
    
    initdensity_flag=0;
    
    lnpdf_init=zeros(size(X,1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now begin tempered smc steps towards next k%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

continue_flag=1;

gamma=0;

while continue_flag
    
    
         
    if initdensity_flag
        
        l_inc=(l+l_nextk)-lnpdf_init;
    else
        l_inc=l_nextk; 
    end
           
    %now find the new_gamma that makes the ESS equal to ESS_bound
    
    %argument in objfun is log(gammaincrement)
      
    log_gamma_inc=fzero('f_Objfun_newgamma',log(.05),[],lnw,l_inc,smcsettings.ESS_bound);
    
    gammanew=min(1,gamma+exp(log_gamma_inc));
                
    incrementalw=(gammanew-gamma)*l_inc;
    
    j=find(isnan(incrementalw) | imag(incrementalw)~=0);
    
    incrementalw(j)=-inf;
            
    lnw=lnw+incrementalw;
    
    %compute normalized weights and compute ESS
    W=exp(lnw-max(lnw));
    
    W=W/sum(W);
    
    ESS=1/sum(W.^2);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %resample and move if ESS<ESS_bound%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    if ESS<smcsettings.ESS_bound
        
        count_RM=count_RM+1;
        
        if smcsettings.verbose
        
            disp(['Resample at k=' num2str(k_current+1) ', initdensity_flag=' num2str(initdensity_flag) ', gamma=' num2str(gammanew)]);
                        
        end
        
        %%%%%%%%%%
        %Resample%
        %%%%%%%%%%
        [X,lnw,l,l_nextk,lnpdf_init]=f_ResampleSet(X,lnw,l,l_nextk,lnpdf_init);
                   
        %%%%%%
        %Move%
        %%%%%%
        
        CumAcceptRate=0;
        
        counter=0;
        
        while (CumAcceptRate<smcsettings.CumAcceptRateBound)
            
            counter=counter+1;
            
            [X,l,l_nextk,lnpdf_init,AcceptRate,filtersettings,smcsettings]=...
                f_MoveSet(X,lnpdf_init,l,l_nextk,k_current,dataX,dataY,gamma,filtersettings,smcsettings,initdensity_flag);
            
            CumAcceptRate=CumAcceptRate+AcceptRate;
            
            if smcsettings.verbose
                
               disp(['Cumulate acceptance rate is ' num2str(CumAcceptRate ) ' in iteration ' num2str(counter)]);
                
            end
                                     
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %If Acceptance Rate falls below a pre-set value, increase the number
        %of particles and re-set source density to the Initialization Density
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if AcceptRate<smcsettings.ChangeMBound
            
            %reinitialize particle poulation at a fitted normal
            initdensity_flag=1;
                        
            smcsettings.Initialization_param=feval(smcsettings.InitializationFit,X);
            
            Nparam=size(X,1);
            
            X=feval(smcsettings.InitializationSim,smcsettings.Initialization_param,Nparam);
            
            lnpdf_init=feval(smcsettings.InitializationLogl,X,smcsettings.Initialization_param);
            
            %double the number of particles
            filtersettings.Nparticles=filtersettings.Nparticles*2;
            
            if smcsettings.verbose
                
                disp(['Nparticles doubled to ' num2str( filtersettings.Nparticles)]);
                
            end
            
            %compute necessary loglikelihood values at X
            l=zeros(Nparam,1);
            
            if k_current>0
                
                for j=1:k_current
                    
                    [l_inc,filtersettings]=feval(filtersettings.filterfun,X,dataX,dataY,filtersettings);
                    
                    l=l+l_inc;
                end
                
            end
            
            [l_nextk,filtersettings]=feval(filtersettings.filterfun,X,dataX,dataY,filtersettings);
            
            lnw = zeros(Nparam,1);
            
            gammanew=0;
                                   
        end
        
                    
    end
    
    
    gamma=gammanew;
        
    continue_flag=(gamma<1);
    
end
  
l=l+l_nextk;

if count_RM==0
    AcceptRate=NaN;
end
    
    





