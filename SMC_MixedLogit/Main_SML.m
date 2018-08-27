%run SML

clear all

addpath '.\mlogit'
randn('state',0);
rand('state',0);

%prepare the data
load data_fish
nprods=4;  % # of alternatives
ni=size(X,1)/nprods;
nbeta=size(X,2);
dataX=reshape(X,[nprods,ni,nbeta]);
dataY=find(Y)-(0:nprods:nprods*(ni-1))';
clear X Y

myhalton = haltonset(2,'Skip',1e3,'Leap',1e2);
myepsilon=norminv(net(myhalton,10000*100));

myobj=@(X)-1*f_Get_logl_MixedLogit_fixv(X,dataX,dataY,myepsilon,10000);%obj for minimization



mystart(:,3:4)=log(mystart(:,3:4)); %reparametrize the sd
myoptions = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','MaxFunctionEvaluations',10000,'MaxIterations',1000);



for i=1:50
    
    i
 tic
  [output_mle.solution(i,:),output_mle.fval(i),output_mle.exitflag(i)]=fminunc(myobj,mystart(i,:),myoptions);
 output_mle.runtime(i)=toc;
end

output_mle.solution(:,3:4)=exp(output_mle.solution(:,3:4));%reparametrize the sd
output_mle.fval= -output_mle.fval; %convert it back to likelihood
save result_mlogit_mle output_mle 
