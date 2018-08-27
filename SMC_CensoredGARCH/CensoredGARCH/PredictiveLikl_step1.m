function loglik_y=PredictiveLikl_step1(Pnew,Pold,mu,beta,Bounds,h,lnP_intrinsic)

%how much rounding we allow around the bounds
Margin=.001;

%find Bounds if any, Bounds expressed in simple returns
Ret=Pnew/Pold-1;

jbound=find(and(abs(Ret)>Bounds-Margin,abs(Ret)<Bounds+Margin));

if ~isempty(jbound)
   
    %case of bound violation
        
    if Ret<0
        
        %active bound is lower bound and is i nterms of log returns the
        %following
        LB=log(1-Bounds(jbound));
        
        loglik_y=log(GEDSTD_cdf((log(Pold)-lnP_intrinsic+LB-mu)./sqrt(h),beta));
                
    else
        %active bound is upper bound and is i nterms of log returns the
        %following
        UB=log(1+Bounds(jbound));
        
        loglik_y=log(1-GEDSTD_cdf((log(Pold)-lnP_intrinsic+UB-mu)./sqrt(h),beta));
          
    end
        
else
    
    %case of no bound violation
            
    loglik_y=GEDSTD_lnpdf((log(Pnew)-lnP_intrinsic-mu)./sqrt(h),beta)-log(h)/2;
      
end

    


          