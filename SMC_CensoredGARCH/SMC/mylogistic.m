function f=mylogistic(x,a,b)

Nparam=size(x,1);


f=a+(b-a)./(1+exp(-x));