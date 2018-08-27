function X=GEDSTD_invcdf(p,beta)
%returns the inverse cdf of standardized ged distribution

%scale parameter ensuring a variance of 1
alpha=sqrt(gamma(1./beta)./gamma(3./beta));

X  = GEDinv(p,alpha,beta);
