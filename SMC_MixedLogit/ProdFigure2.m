clear all

load result_mlogit_autop;

K=size(output(1).Xmean,1);

Nsim=length(output);

logl=zeros(Nsim,K);

Nparticles=zeros(Nsim,K);

runtime=zeros(Nsim,K);

for i=1:Nsim
    
    Nparticles(i,:)=output(i).Nparticles;
    
    runtime(i,:)=cumsum(output(i).runtime);
    
   logl(i,:)=output(i).loglik;
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

xlabel('m');

title('Mean number of latent state particles');

subplot(3,1,3);

plot(mean(runtime));

xlim([1 K]);

xlabel('m');

title('Mean runtime in sec');

print -dpdf Figure_MLogit.pdf


