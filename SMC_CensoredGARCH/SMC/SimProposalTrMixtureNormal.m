function [X_proposal,lnprop_X,lnprop_Xprop]=SimProposalTrMixtureNormal(X,expindex,logisticindex,...
    logisticbound_l,logisticbound_u,mixturedim)
%simulate from mv normal fitted on X allows transformed marginals
%(exp,logistic)
%assume even-weighted sample

warning off;

[Nparam,K]=size(X);

%get transformed X

X_tr=X;


if ~isempty(expindex)
   
    X_tr(:,expindex)=log(X_tr(:,expindex));
    
end

if ~isempty(logisticindex)
    
    logisticbound_l=repmat(logisticbound_l,Nparam,1);

    logisticbound_u=repmat(logisticbound_u,Nparam,1);
    
    X_tr(:,logisticindex)=mylogisticinv(X_tr(:,logisticindex),logisticbound_l,logisticbound_u);
    
end

%fit mv mixture normal on X_tr and sim from it

options = statset('Display','off');

try 

     %fit gaussian mixture 
     	obj = gmdistribution.fit(X_tr,mixturedim,'Regularize',1e-6,'Options',options);

catch
 
     	obj = gmdistribution.fit(X_tr,mixturedim,'SharedCov',true,'Options',options);
end

X_tr_proposal=random(obj,Nparam);

lng_prop=log(pdf(obj,X_tr_proposal));

lng=log(pdf(obj,X_tr));

%compute jacobian
logJacobian_prop=zeros(Nparam,1);

logJacobian=zeros(Nparam,1);

if ~isempty(expindex)

logJacobian_prop=logJacobian_prop+...
    sum(-X_tr_proposal(:,expindex),2);


logJacobian=logJacobian+...
    sum(-X_tr(:,expindex),2);

end


if ~isempty(logisticindex)

logJacobian_prop=logJacobian_prop+...
    sum(-log(mylogisticder(X_tr_proposal(:,logisticindex),logisticbound_l,logisticbound_u)),2);


logJacobian=logJacobian+...
    sum(-log(mylogisticder(X_tr(:,logisticindex),logisticbound_l,logisticbound_u)),2);

end

lnprop_X=lng+logJacobian;

lnprop_Xprop=lng_prop+logJacobian_prop;

%transform proposal

X_proposal=X_tr_proposal;

if ~isempty(expindex)
   
    X_proposal(:,expindex)=exp(X_tr_proposal(:,expindex));
    
end

if ~isempty(logisticindex)
   
    X_proposal(:,logisticindex)=mylogistic(X_tr_proposal(:,logisticindex),logisticbound_l,logisticbound_u);
    
end
        


