function initial_sample=f_PriorSim_MixedLogit(prior_param,Nparam)
% alpha delta sigmaV
%
% use normal initialization alpha and truncated normal for delta and sigmaV
% 2/28/16 for mixed logit
% 5/13/16 non-negativity 4:6

% Trunc normal
initial_sample=repmat(prior_param.mu,Nparam*5,1)+...
              repmat(prior_param.sig,Nparam*5,1).*randn(Nparam*5,6);

      
          
rowindex=find(initial_sample(:,4)>=0 & initial_sample(:,5)>=0 & initial_sample(:,6)>=0);

initial_sample=initial_sample(rowindex(1:Nparam,:),:);

% Normal
%  initial_sample=repmat(prior_param.mu,Nparam,1)+...
%                repmat(prior_param.sig,Nparam,1).*randn(Nparam,6);


