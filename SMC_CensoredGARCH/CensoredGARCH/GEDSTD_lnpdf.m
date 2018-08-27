function F=GEDSTD_lnpdf(X,beta)
%returns the log pdf of standardized ged distribution

%scale parameter ensuring a variance of 1
alpha=sqrt(gamma(1./beta)./gamma(3./beta));

F=log(beta)-log(2)-log(alpha)-gammaln(1./beta)-(abs(X)./alpha).^beta;