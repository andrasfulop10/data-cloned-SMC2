% %add folder containing the model specific filtering files and SMC code
 clear all
 
 addpath '.\CensoredGARCH';
 
 addpath '.\SMC';
 
 randn('state',0);
 rand('state',0);
 
%%%%%%%%
%run smc%
%%%%%%%%%

clear all;

load ChineseData;

%loop over firms
for ifirms=1:100

    
filtersettings.filterfun='Get_logl';

%number of periods that is used to simulate variance state variable
%initially to get to stationary distribution, which is then used to initialize
%the particle filter.
filtersettings.Ti=50;

%number of replicates per particle when there is a data gap between two
%observations
filtersettings.Nrepgap=5;

%this parameter defines when does the code is parallelized in the clone
%dimensions. The reason is to keep the memory blueprint of the algorithm
%under control
%if Nparticles<filtersettings.cloneparallelprocessNparticles, the code is
%parallelized in the clone direction ow, it is not.
%this may lead the algorithm to shut down this parallelization as
%Nparticles increases at higher clone numbers
filtersettings.cloneparallelprocessNparticles=1600;

%number of individual loglikelihoods that are dropped at the beginning of
%the data sample (objective is to make teh objective function more robust
%to the initialization of the hidden state variable)
filtersettings.burnin=10;

%filtersettings.gpu_switch=1 means that the code uses a gpu to speed up the
%likelihood vealuations. This necessitates Matlab parallel toolbox
%installed and an Nvidia cuda-capable gpu.
%filtersettings.gpu_switch=0 means only cpu is used
filtersettings.gpu_switch=1;

%possible price limits
filtersettings.Bounds=[0.05 0.1];

%number of parameter particles 
Nparam=1000;

%smc stop 
%if smcsettings.Kend , (maximum number of clones) is achieved
smcsettings.Kend=10;

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
smcsettings.mixturedim=4;

%acceptance rate below which number of state particles is doubled
smcsettings.ChangeMBound=.1;

%maximum number of state particles. this is needed to keep the memory
%blueprint of the algorithm under control
smcsettings.NparticlesUpperBound=6400;

%these two fields only play a role if RW proposal
smcsettings.scale=1;
smcsettings.adaptproposalscale=0;

smcsettings.expindex=[];
smcsettings.logisticindex=[];
smcsettings.logisticbound_l=[];
smcsettings.logisticbound_u=[];

%initialization sampler 
smcsettings.InitializationSim='PriorSim';
smcsettings.InitializationLogl='PriorLogl';
smcsettings.InitializationFit='PriorFit';

%initialization sampler params
Eh=0.3^2*1/250;

alpha1=.7;
alpha2=.1;
theta=1;
mu=0;

beta=1;

alpha0=Eh.*(1-alpha1-alpha2*(1+theta^2));

smcsettings.Initialization_param.mu=InvTransformParam([alpha0 alpha1 alpha2 theta mu beta]);
smcsettings.Initialization_param.V=diag([ones(4,1); .01; 1]);

Data.P=output(ifirms).Price;

%Data.timestamp=(1:length(Data.P))';

Data.timestamp=output(ifirms).timestamp;
 
%initial number of state particles
filtersettings.Nparticles=100;
    [res(ifirms).Xmean,Xsample,res(ifirms).runtime,ESS,res(ifirms).AcceptRate,...
        res(ifirms).Nparticles]...
        =RunSMC(Nparam,filtersettings,smcsettings,Data);
    
    [ph1,ph2,res(ifirms).FiltMeans]=...
        Get_logl(res(ifirms).Xmean(end,:),Data,filtersettings);
    
    save ChineseRes_CensoredGARCH output res Nparam smcsettings filtersettings;   

end


    
    



