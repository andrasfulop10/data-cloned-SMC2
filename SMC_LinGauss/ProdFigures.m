clear all

addpath '.\LinGaussFilter'

load simres_LinGauss;

K=size(output(1).Xmean,1);

Nsim=length(output);

logl=zeros(Nsim,K);

Nparticles=zeros(Nsim,K);

runtime=zeros(Nsim,K);

for i=1:Nsim
    
    Nparticles(i,:)=output(i).Nparticles;
    
    runtime(i,:)=cumsum(output(i).runtime);
    
    for l=1:K
    
        X=output(i).Xmean(l,:);
        
        [logl(i,l),filtersettings]=KF_logl(X,Y,filtersettings);
            
    end
end

subplot(3,1,1);

for i=1:Nsim
   
    plot(1:K,logl(i,:),'bl');
    
    hold on;
    
end

hold off;

xlim([1 K]);

xlabel('m');

title('LogLikelihood');

subplot(3,1,2);

plot(mean(Nparticles));

xlim([1 K]);

ylim([0 800]);

yticks(0:200:800);

xlabel('m');

title('Mean number of latent state particles');

subplot(3,1,3);

plot(mean(runtime));

xlim([1 K]);

xlabel('m');

ylim([0 3000]);

yticks(0:1000:3000);

title('Mean runtime in sec');

print -dpdf Figure_LinGauss.pdf

Kvec=[1 5 10 20];

for i=1:length(Kvec)
   
    c=[num2str(Kvec(i)) ' & ' num2str(mean(logl(:,Kvec(i))))  ' & ' num2str(std(logl(:,Kvec(i)))) ...
        ' & '  num2str(mean(runtime(:,Kvec(i))))     ' \\ '];
        
    disp(c);
       
end


