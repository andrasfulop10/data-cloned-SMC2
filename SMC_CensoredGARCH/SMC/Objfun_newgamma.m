function f=Objfun_newgamma(log_gamma_inc,lnw,l_inc,ESS_bound)


incrementalw=exp(log_gamma_inc)*l_inc;

j=find(isnan(incrementalw) | imag(incrementalw)~=0);

incrementalw(j)=-inf;

lnw=lnw+incrementalw;

W=exp(lnw-max(lnw));

ESS=sum(W).^2./sum(W.^2);

f=ESS-(ESS_bound-1);


