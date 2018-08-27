function l=KF_objfun(X,Y,filtersettings)

l=KF_logl(X',Y,filtersettings);

l=-l;