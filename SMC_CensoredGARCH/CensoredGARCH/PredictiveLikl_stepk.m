function  [loglik_y,lnP_intrinsic_beforet,h_beforet]=...
            PredictiveLikl_stepk(Pnew,alpha0,alpha1,alpha2,theta,mu,beta,...
                Bounds,h,lnP_intrinsic,Nrepgap,timegap,gpu_switch)
 
[Nparticles,Nparam]=size(h);            
     
%allocate output matrices
if gpu_switch
    
     loglik_y=gpuArray.zeros(Nparticles*Nrepgap,Nparam);
     lnP_intrinsic_beforet=gpuArray.zeros(Nparticles*Nrepgap,Nparam);
     h_beforet=gpuArray.zeros(Nparticles*Nrepgap,Nparam);
          
 else
     loglik_y=zeros(Nparticles*Nrepgap,Nparam);
     lnP_intrinsic_beforet=zeros(Nparticles*Nrepgap,Nparam);
     h_beforet=zeros(Nparticles*Nrepgap,Nparam);
          
end
 
 %run Nrepgap loops to fill up output matrices by forward simulation and
 %likelihood evaluation
 for isim=1:Nrepgap
    
     %%%%%%%%%%%%%%%%%%%%%%%%%
     %forward simulation loop%
     %%%%%%%%%%%%%%%%%%%%%%%%%
     
     %initialize
     hcurrent=h;
     
     if gpu_switch
        cumret=gpuArray.zeros(Nparticles,Nparam);
     else
        cumret=zeros(Nparticles,Nparam);
     end
     
     %roll forward
     for t=1:timegap-1
         
         if gpu_switch
             u=gpuArray.rand(Nparticles,Nparam);
         else
             u=rand(Nparticles,Nparam);
         end
         
         epsilon=GEDSTD_invcdf(u,beta);
              
         cumret=mu+sqrt(hcurrent).*epsilon;
         
         hcurrent=alpha0+alpha1.*hcurrent+alpha2.*hcurrent.*(epsilon-theta).^2;
          
     end
           
     colindex=(1:Nparticles)+(isim-1)*Nparticles;
     
     lnP_intrinsic_beforet(colindex,:)=lnP_intrinsic+cumret;
     
     h_beforet(colindex,:)=hcurrent;
         
     loglik_y(colindex,:)=GEDSTD_lnpdf...
         ((log(Pnew)-lnP_intrinsic_beforet(colindex,:)-mu)./sqrt(h_beforet(colindex,:)),beta)...
         -log(h_beforet(colindex,:))/2;
               
 end
 
 