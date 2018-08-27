function [l,E_alpha,sig2_alpha]=KF(mu,phi,sigma2_eta,sigma2_epsilon,Y,E_alpha,sig2_alpha)

T=length(Y);

Nparam=size(mu,1);

l=zeros(Nparam,T);

for t=1:T
    
    %%%%%%%%%%%%%%
    %update state%
    %%%%%%%%%%%%%%
    
    E_alpha=mu+phi.*E_alpha;
    
    sig2_alpha=phi.^2.*sig2_alpha+sigma2_eta;
    
    %%%%%%%%%%%%%%%%%%%%%
    %digest observation%
    %%%%%%%%%%%%%%%%%%%%%
    
    %predicted moments of observation
    Ey=E_alpha;
    
    sig2y=sig2_alpha+sigma2_epsilon;
    
    %evaluate log likelihood
    l(:,t)=-(Y(t)-Ey).^2./(2*sig2y)-log(sig2y)/2-log(2*pi)/2;
    
    %update state
    
    E_alpha=E_alpha+sig2_alpha./sig2y.*(Y(t)-Ey);
    
    sig2_alpha=sig2_alpha.*(1-sig2_alpha./sig2y);
    
end



