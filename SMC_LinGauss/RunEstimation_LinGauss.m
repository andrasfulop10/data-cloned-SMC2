%add folder containing the model specific filtering files
clear all

addpath '.\LinGaussFilter'

randn('state',0);
rand('state',0);

%%%%%%%%%%%%%%%%%%%%%
%simulate from model%
%%%%%%%%%%%%%%%%%%%%%

%real parameters same as in flurry and shepard
mu=0.5;
sigma_epsilon=sqrt(1);
phi=.825;
sigma_eta=sqrt(.75);

T=500;

%simulate data

alpha=zeros(T,1);
Y=zeros(T,1);

%simulate from stationary distribution of hidden variable
alpha_current=sigma_eta/sqrt(1-phi^2)*randn;

%run loop to simulate from model

for t=1:T
   
    alpha_current=phi*alpha_current+sigma_eta*randn;
    
    alpha(t)=alpha_current;
    
    Y(t)=mu+alpha_current+sigma_epsilon*randn;
    
end

%%%%%%%%
%run smc%
%%%%%%%%%

%filtersettings.filterfun is the name of the Matlab function used for
%likelihood evaluation
%filtersettings.filterfun='KF_logl';

filtersettings.filterfun='PF_adapted_logl';

%filtersettings.gpu_switch=1 means that the code uses a gpu to speed up the
%likelihood vealuations. This necessitates Matlab parallel toolbox
%installed and an Nvidia cuda-capable gpu.
%filtersettings.gpu_switch=0 means only cpu is used
filtersettings.gpu_switch=1;

% number of fixed parameter sets in SMC
Nparam=500;

%filtersettings.randomnumbers.CRN=1: common random numbers used in
%likelihood evaluations across parameter sets. 
%filtersettings.randomnumbers.CRN=0: likelihood evaluations across
%parameter sets should use independent random numbers. This should be the
%case to be in the pseudo-marginal setting
filtersettings.randomnumbers.CRN=0;

%smc stop 
%if smcsettings.Kend , (maximum number of clones) is achieved
smcsettings.Kend=20;

%amount of info SMC algorithm displays
smcsettings.verbose=1;

%smcsettings.ESS_bound: when is the resample-move step triggered
smcsettings.ESS_bound=Nparam/2;

%smcsettings.CumAcceptRateBound: cumulative acceptance rate when the move
%steps are stopped
smcsettings.CumAcceptRateBound=1;

%proposal used in move step: 
%smcsettings.proposal='RW' random walk
%smcsettings.proposal='MixtureNormal' Independent Mixture of normals
smcsettings.proposal='MixtureNormal';

%number of mixture components in Independent Mixture proposal
smcsettings.mixturedim=1;

%acceptance rate below which number of state particles is doubled
smcsettings.ChangeMBound=.2;

%these two fields only play a role if random walk proposal
smcsettings.scale=1;
smcsettings.adaptproposalscale=0;

%variable transformation in proposal to automatically satisfy positiity
%(exp) or boundedness (logistic)
smcsettings.expindex=[];
smcsettings.logisticindex=[2];
smcsettings.logisticbound_l=[-1];
smcsettings.logisticbound_u=[1];

%initialization sampler 
smcsettings.InitializationSim='PriorSim';
smcsettings.InitializationLogl='PriorLogl';
smcsettings.InitializationFit='PriorFit';
smcsettings.Initialization_param.mu=[.25 0 log(1.5) .475];
smcsettings.Initialization_param.V=eye(4);

%number of independent runs in SMC
Nsim=50;

for isim=1:Nsim
    
    %initial number of state particles
    filtersettings.Nparticles=50;
    
    [output(isim).Xmean,Xsample,output(isim).runtime,ESS,output(isim).AcceptRate,...
        output(isim).Nparticles]...
        =RunSMC(Nparam,filtersettings,smcsettings,Y);

 save simres_LinGauss output Y Nparam smcsettings filtersettings;   

end


    
    



