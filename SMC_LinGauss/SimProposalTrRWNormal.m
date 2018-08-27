function [X_proposal,lnprop_X,lnprop_Xprop]=SimProposalTrRWNormal(X,expindex,logisticindex,...
    logisticbound_l,logisticbound_u,scale)
%simulate from mv normal fitted on X allows transformed marginals
%(exp,logistic)
%assume even-weighted sample

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

%fit mv normal on X_tr and sim from it

V_proposal=cov(X_tr);

%draw from epsilon using eigen decomposition

[U,Lambda]=eig(V_proposal); %U*Lamda*U'=V_proposal

ld=diag(Lambda);

j=find(Lambda<max(ld)*1e-6);

Lambda(j)=0;

Keff=sum(diag(Lambda)>0);

Sig=U*Lambda.^.5;

z=(Sig*randn(K,Nparam))';

epsilon=2.38/sqrt(Keff)*scale*z;

X_tr_proposal=X_tr+epsilon;

%lng and lng_prop cancel out as RW proposal
lng_prop=zeros(Nparam,1);

lng=zeros(Nparam,1);

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
        


