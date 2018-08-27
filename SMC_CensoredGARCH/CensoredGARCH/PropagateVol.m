function [h,lnP_intrinsic]=PropagateVol(h,lnP_intrinsic,Pnew,Pold,alpha0,alpha1,alpha2,theta,mu,beta,...
    Bounds,gpu_switch)

[Nparticles,Nparam]=size(h);

%how much rounding we allow around the bounds
Margin=.001;

%find Bounds if any, Bounds expressed in simple returns
Ret=Pnew/Pold-1;

jbound=find(and(abs(Ret)>Bounds-Margin,abs(Ret)<Bounds+Margin));

if ~isempty(jbound)
    
    %case of bound violation
    if gpu_switch
        randu=gpuArray.rand(Nparticles,Nparam);
    else
        randu=rand(Nparticles,Nparam);
    end
    
    if Ret<0
        
        %active bound is lower bound and is i nterms of log returns the
        %following
        LB=log(1-Bounds(jbound));
        
        L_PR=GEDSTD_cdf((log(Pold)-lnP_intrinsic+LB-mu)./sqrt(h),beta);
        
        epsilon=GEDSTD_invcdf(randu.*L_PR,beta);
        
    else
        %active bound is upper bound and is i nterms of log returns the
        %following
        UB=log(1+Bounds(jbound));
        
        U_PR=GEDSTD_cdf((log(Pold)-lnP_intrinsic+UB-mu)./sqrt(h),beta);
        
        epsilon=GEDSTD_invcdf(U_PR+(1-U_PR).*randu,beta);
                
    end
    
else
    
    %case of no bound violation
    epsilon=(log(Pnew)-lnP_intrinsic-mu)./sqrt(h);
    
end

%equation for intrinsic price
lnP_intrinsic=lnP_intrinsic+mu+sqrt(h).*epsilon;

%GARCH equation
h=alpha0+alpha1.*h+alpha2.*h.*(epsilon-theta).^2;


            
            