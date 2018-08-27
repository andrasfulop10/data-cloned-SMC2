clear all;

randn('state',0);

rand('state',0);

%real parameters same as in flurry and shepard
mu=0.5;
sigma_epsilon=sqrt(1);
phi=.825;
sigma_eta=sqrt(.75);

T=1000;

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

%run filters
Nparam=100;

%parameterization X= [mu phi log(sigma_eta) log(sigma_epsilon)]

X=repmat([mu phi log(sigma_eta) log(sigma_epsilon)],Nparam,1);

filtersettings=[];

%run Kalman Filter
[l,filtersettings]=KF_logl(X,Y,filtersettings);

%run adapted PF
filtersettings.gpu_switch=0;

filtersettings.Nparticles=100;

filtersettings.randomnumbers.CRN=1;

filtersettings.randomnumbers.epsilon=randn(filtersettings.Nparticles,T+1);
filtersettings.randomnumbers.u=rand(filtersettings.Nparticles,T);

[l_PF,filtersettings]=PF_adapted_logl(X,Y,filtersettings);













