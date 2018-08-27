function f=mvnlnpdf2(x,mu,sigma)
 %compute the log-density of multivariate normal
 %mu must be row vector
 f=log(mvnpdf(x,mu,sigma));
 