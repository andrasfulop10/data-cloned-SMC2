function [X,l,l_nextk,lnpdf_init,AcceptRate,filtersettings,smcsettings]=...
    f_MoveSet(X,lnpdf_init,l,l_nextk,k_current,dataX,dataY,gamma,filtersettings,smcsettings,initdensity_flag)
           
%function assumes an equal-weighted param set

Nparam=size(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%
%simulate from proposal%
%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(smcsettings.proposal,'RW')
    
    %normal random walk proposal
    % use exponential transformations for smcsettings.expindex params to
    % keep them positive
    % use logistic transformations for smcsettings.logisticindex params to
    % keep them in a bounded domain
    
    [X_proposal,logprop_orig,logprop_proposal]=...
        f_SimProposalTrRWNormal(X,smcsettings.expindex,smcsettings.logisticindex,...
        smcsettings.logisticbound_l,smcsettings.logisticbound_u,smcsettings.scale);
    
    
elseif strcmp(smcsettings.proposal,'MixtureNormal')
    %independent normal mixture proposal
    % use exponential transformations for smcsettings.expindex params to
    % keep them positive
    % use logistic transformations for smcsettings.logisticndex params to
    % keep them in a bounded domain
    %
    %NOTE: For now the code does not deal with linear or nonlinear constraints between the parameters. 
    %However these could be easily included by extending
    %SimProposalTrMixtureNormal such that it only proposes within the
    %allowed domain. Evaluating the probability of the constrained set is
    %not necessary, as it cancels out in the acceptance rates for
    %independent proposals (not for the RW sampler though!!!)
    
    [X_proposal,logprop_orig,logprop_proposal]=...
        f_SimProposalTrMixtureNormal(X,smcsettings.expindex,smcsettings.logisticindex,...
        smcsettings.logisticbound_l,smcsettings.logisticbound_u,smcsettings.mixturedim);
else
    
    error('Unknown Proposal Distribution in Move Step');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute initialization sampler at proposal%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initdensity_flag

    lnpdf_init_proposal=feval(smcsettings.InitializationLogl,X_proposal,smcsettings.Initialization_param);
else
    lnpdf_init_proposal=zeros(Nparam,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimated data-cloned loglikelihood at proposal            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l_proposal=zeros(Nparam,1);

if k_current>0
    
    for j=1:k_current
        
        [l_inc,filtersettings]=feval(filtersettings.filterfun,X_proposal,dataX,dataY,filtersettings);
        
        l_proposal=l_proposal+l_inc;
    end
    
 end
    
 [l_nextk_proposal,filtersettings]=feval(filtersettings.filterfun,X_proposal,dataX,dataY,filtersettings);
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute acceptance weights%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initdensity_flag
    
    logtargetn_proposal=(l_nextk_proposal+l_proposal)*gamma+(1-gamma)*lnpdf_init_proposal;
    
    logtargetn_orig=(l_nextk+l)*gamma+(1-gamma)*lnpdf_init;
    
else
    
    logtargetn_proposal=l_proposal+l_nextk_proposal*gamma;
    
    logtargetn_orig=l+l_nextk*gamma;
    
end

j=find(imag(logtargetn_proposal)~=0 | ~isfinite(logtargetn_proposal));

logtargetn_proposal(j)=-inf;

lnalpha=logtargetn_proposal-logtargetn_orig+logprop_orig-logprop_proposal;

logu=log(rand(Nparam,1));

accept=find(logu<lnalpha);

AcceptRate=length(accept)/Nparam;

if smcsettings.verbose
    
    disp(['Acceptance rate in MH step: ' num2str(AcceptRate) 'at scale ' num2str(smcsettings.scale)]);
end

%adapt scale of proposal if necessary
if and(strcmp(smcsettings.proposal,'RW'),smcsettings.adaptproposalscale)
    
    if AcceptRate<.2
        smcsettings.scale=smcsettings.scale*.75;
    elseif AcceptRate>.4
        smcsettings.scale=smcsettings.scale/.75;
    end
    
end

%implement move

X(accept,:)=X_proposal(accept,:);

l(accept)=l_proposal(accept);

l_nextk(accept)=l_nextk_proposal(accept);

lnpdf_init(accept)=lnpdf_init(accept);


