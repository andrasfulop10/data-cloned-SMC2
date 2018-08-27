function Xtransform=TransformParam(X)

%alpha0
alpha0=exp(X(:,1));

%theta
theta=X(:,4);

%alpha1,alpha2
h1=mylogistic(X(:,2),0,1);

h2=mylogistic(X(:,3),0,1);

alpha1=h1.*h2;

alpha2=(h1-alpha1)./(1+theta.^2);

mu=X(:,5);

beta=mylogistic(X(:,6),0,2);

Xtransform=[alpha0 alpha1 alpha2 theta mu beta];
