function F=GEDSTD_cdf(X,beta)
%returns the cdf of standardized ged distribution

%scale parameter ensuring a variance of 1
alpha=sqrt(gamma(1./beta)./gamma(3./beta));

F=GEDcdf( X, alpha, beta );