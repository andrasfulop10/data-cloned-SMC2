function [X,lnw,l,l_nextk,lnpdf_prior]=f_ResampleSet(X,lnw,l,l_nextk,lnpdf_prior)
%this function resamples the fixed parameter set and the quantities needed
%for the algorithm

Nparam=size(X,1);

cumw = [0; cumsum(exp(lnw - max(lnw)))];
cumw = cumw/cumw(end);

%use stratified resampling
[PH, bin] = histc(((0:Nparam-1)+rand(1,Nparam))/Nparam,cumw);
X = X(bin,:);
l = l(bin);
l_nextk=l_nextk(bin);

lnpdf_prior = lnpdf_prior(bin);

lnw = zeros(Nparam, 1);
