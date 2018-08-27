function [Xmean,runtime]=f_main_test_SMC2(start,Nparam,Nparticles,Nclone,dataX,dataY,index_logn,opt_fixv,opt_smc2)
%4/29/16 compared with Main_Test_SMC, change the following:
%        filtersettings.Nparticles=Nparticles
%        smcsettings.Kend=Nclone   
%        add folder, random seed setting outside
%        data loading outside
%        epsilon generated outside
%        Nparam=6000 outside
%        smcsettings.mixturedim=2
%5/13/16 smcsettings.expindex=4:6
%        smcsettings.mixturedim=1
%        add index_logn, if =1
%        set smcsettings.InitializationSim/Logl to option 2
%        set prior mu to log(mu)
%        not change prior sigma
%6/12/16 add input arg: sobol opt_fixv opt_smc2
        
%add folder containing the model specific filtering files

% Moved outside to Main_New

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data: instead of simulation, load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Moved outside to Main_New

%%%%%%%%
%run smc%
%%%%%%%%%

%filtersettings.filterfun='KF_logl'; 'PF_adapted_logl';

if opt_fixv
    filtersettings.filterfun='f_Get_logl_MixedLogit_fixv'; 
else
    filtersettings.filterfun='f_Get_logl_MixedLogit_randv'; 
end


filtersettings.gpu_switch=0;

filtersettings.Nparticles=Nparticles;

% Nparam=Nparam_vec(i_N);


smcsettings.SMC2=opt_smc2;

if smcsettings.SMC2==1
    
    filtersettings.randomnumbers.CRN=0;
else
    filtersettings.randomnumbers.CRN=1;
%     filtersettings.randomnumbers.epsilon=randn(filtersettings.Nparticles,T+1);
%     filtersettings.randomnumbers.u=rand(filtersettings.Nparticles,T);
    
end

smcsettings.Kend=Nclone;                   
smcsettings.verbose=1;
smcsettings.ESS_bound=Nparam/2;
smcsettings.ESS_bound_Resample=Nparam/1;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smcsettings.CumAcceptRateBound=1;
smcsettings.ChangeMBound=0.2;
smcsettings.proposal='MixtureNormal';

%this field only plays a role in mixutre proposal
smcsettings.mixturedim=7;

%these two fields only play a role if RW proposal
smcsettings.scale=1;
smcsettings.adaptproposalscale=0;

smcsettings.expindex=3:4;              % Index of logNormal for proposal
smcsettings.logisticindex=[];          % Index of uniform for proposal
smcsettings.logisticbound_l=[];
smcsettings.logisticbound_u=[];

%initialization sampler 
smcsettings.InitializationFit='PriorFit';
smcsettings.Initialization_param.mu=start;
smcsettings.Initialization_param.sig=3*eye(4);

%if index_logn==0
% smcsettings.InitializationSim='f_PriorSim_MixedLogit';
% smcsettings.InitializationLogl='f_PriorLogl_MixedLogit';
%end

if  index_logn==1
    smcsettings.Initialization_param.mu(:,3:4)=log(smcsettings.Initialization_param.mu(:,3:4));
    smcsettings.InitializationSim='f_PriorSim_MixedLogit2';
    smcsettings.InitializationLogl='f_PriorLogl_MixedLogit2';
end
   




   
[Xmean,Xsample,runtime,ESS,AcceptRate]=f_RunSMC(Nparam,filtersettings,smcsettings,dataX,dataY);

