% %add folder containing the model specific filtering files adn SMC code
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

for ifirms=1:100

filtersettings.filterfun='Get_logl_GARCH';

filtersettings.cloneparallelprocessNparticles=400;

filtersettings.burnin=10;

filtersettings.gpu_switch=0;

Nparam=1000;

%smc stop 
%either if we cannot reject the null that increase in
%loglikelihood is greater than  smcsettings.stopcriterion
%or if smcsettings.Kend is achieved
smcsettings.Kend=10;
%smcsettings.stopcriterion=.1;

smcsettings.verbose=1;
smcsettings.ESS_bound=Nparam/2;
smcsettings.CumAcceptRateBound=1;

smcsettings.proposal='MixtureNormal';

%this field only plays a role in mixture proposal
smcsettings.mixturedim=4;

%acceptance rate below which number of state particles is doubled
smcsettings.ChangeMBound=0;

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
Data.timestamp=output(ifirms).timestamp;
    
%initial number of state particles
filtersettings.Nparticles=1;

    [res(ifirms).Xmean,Xsample,res(ifirms).runtime,ESS,res(ifirms).AcceptRate,...
        res(ifirms).Nparticles]...
        =RunSMC(Nparam,filtersettings,smcsettings,Data);

    
    [ph1,ph2,res(ifirms).FiltMeans]=...
        Get_logl_GARCH(res(ifirms).Xmean(end,:),Data,filtersettings);
        
    save ChineseRes_GARCH output res Nparam smcsettings filtersettings;   

end


    
    



