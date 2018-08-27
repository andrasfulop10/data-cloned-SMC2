function [h,lnP_intrinsic]=PropagateVol_stepk(Pnew,lnP_intrinsic_beforet,h_beforet,...
    alpha0,alpha1,alpha2,theta,mu,gpu_switch)

[Nparticles,Nparam]=size(h_beforet);
    
epsilon=(log(Pnew)-lnP_intrinsic_beforet-mu)./sqrt(h_beforet);
   
lnP_intrinsic=lnP_intrinsic_beforet+mu+sqrt(h_beforet).*epsilon;

%GARCH equation
h=alpha0+alpha1.*h_beforet+alpha2.*h_beforet.*(epsilon-theta).^2;
            
            