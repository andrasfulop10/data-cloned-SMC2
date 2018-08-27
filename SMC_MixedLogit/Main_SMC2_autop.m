%data-cloning smc2
%automatic choice of p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

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


%likelihood based on fixed halton seq.
myhalton = haltonset(2,'Skip',1e3,'Leap',1e2);
myepsilon=norminv(net(myhalton,10000*100));
myobj=@(X)f_Get_logl_MixedLogit_fixv2(X,dataX,dataY,myepsilon,10000)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% algorithm setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameter setup
N_latent_particles=50;  %initial number of latent particles p
Nrep=size(mystart,1);   %number of MC repetitions
Nclone=100;              %max number of clone M
Nparam=500;             %number of parameter particles

%filtersettings.filterfun='f_Get_logl_MixedLogit_randv'; %use standard random number for marginalization

%this is the version where the loglikelihood code has also been
%parallelized in the fixed params
filtersettings.filterfun='f_Get_logl_MixedLogit_randv_parNpar';


%filtersettings.gpu_switch=1 means that the code uses a gpu to speed up the
%likelihood vealuations. This necessitates Matlab parallel toolbox
%installed and an Nvidia cuda-capable gpu.
%filtersettings.gpu_switch=0 means only cpu is used
filtersettings.gpu_switch=1;


%setting initial number of latent particles (p)
filtersettings.Nparticles=N_latent_particles;

%filtersettings.randomnumbers.CRN=1: common random numbers used in
%likelihood evaluations across parameter sets. 
%filtersettings.randomnumbers.CRN=0: likelihood evaluations across
%parameter sets should use independent random numbers. This should be the
%case to be in the pseudo-marginal setting
filtersettings.randomnumbers.CRN=0;


%smc stop 
%if smcsettings.Kend , (maximum number of clones) is achieved
smcsettings.Kend=Nclone; 


%amount of info SMC algorithm displays
smcsettings.verbose=1;


%smcsettings.ESS_bound: when is the resample-move step triggered
smcsettings.ESS_bound=Nparam/2;

%smcsettings.CumAcceptRateBound: cumulative acceptance rate when the move
%steps are stopped
smcsettings.CumAcceptRateBound=1;

%acceptance rate below which number of state particles is doubled.
smcsettings.ChangeMBound=0.2;

%proposal used in move step: 
%smcsettings.proposal='RW' random walk
%smcsettings.proposal='MixtureNormal' Independent Mixture of normals
smcsettings.proposal='MixtureNormal';


%number of mixture components in Independent Mixture proposal
smcsettings.mixturedim=4; %number of mixture

%these two fields only play a role if RW proposal
smcsettings.scale=1;
smcsettings.adaptproposalscale=0;

smcsettings.expindex=3:4;              % Index of logNormal for proposal
smcsettings.logisticindex=[];          % Index of uniform for proposal
smcsettings.logisticbound_l=[];
smcsettings.logisticbound_u=[];

%initialization sampler 
smcsettings.InitializationFit='PriorFit'; 
smcsettings.InitializationSim='f_PriorSim_MixedLogit2';
smcsettings.InitializationLogl='f_PriorLogl_MixedLogit2';

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                

for isim=1:Nrep
    disp(isim);
    %prior mean
    smcsettings.Initialization_param.mu=mystart(isim,:);
    smcsettings.Initialization_param.mu(:,3:4)=log(smcsettings.Initialization_param.mu(:,3:4));%reparametrize the sd
    
    %prior covariance matrix
    smcsettings.Initialization_param.sig=3*eye(4);

    [output(isim).Xmean,Xsample,output(isim).runtime,ESS,output(isim).AcceptRate,...
        output(isim).Nparticles]...
        =f_RunSMC(Nparam,filtersettings,smcsettings,dataX,dataY);

    output(isim).loglik=myobj(output(isim).Xmean);
    
    save result_mlogit_autop output dataX ...
        dataY Xsample Nparam smcsettings ...
        filtersettings output;   

end





