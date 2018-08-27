function x=f_mylogisticinv(f,a,b)

Nparam=size(f,1);


x=-log((b-a)./(f-a)-1);