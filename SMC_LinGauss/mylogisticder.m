function f=mylogisticder(x,a,b)

Nparam=size(x,1);

f=(b-a).*exp(-x)./(1+exp(-x)).^2;