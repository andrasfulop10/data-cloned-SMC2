function ESS=TempFun(gamma_inc,lnw,lndensity_ratio)
%returns ESS(gamma)
incrementalw=bsxfun(@times,exp(gamma_inc),lndensity_ratio);

j=find(isnan(incrementalw) | imag(incrementalw)~=0);

incrementalw(j)=-inf;

lnw=bsxfun(@plus,lnw,incrementalw);

W=exp(bsxfun(@minus,lnw,max(lnw)));

ESS=sum(W).^2./sum(W.^2);


