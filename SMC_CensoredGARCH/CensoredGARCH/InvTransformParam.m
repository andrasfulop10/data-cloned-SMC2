function X=InvTransformParam(Xtransform)

h0=log(Xtransform(:,1));

theta=Xtransform(:,4);

alpha1=Xtransform(:,2);

alpha2=Xtransform(:,3);

h1=mylogisticinv(alpha1+alpha2.*(1+theta.^2),0,1);

h2=mylogisticinv(alpha1./(alpha1+alpha2.*(1+theta.^2)),0,1);

mu=Xtransform(:,5);

h6=mylogisticinv(Xtransform(:,6),0,2);

X=[h0 h1 h2 theta mu h6];