function f=mvnlnpdf(x,mu,sigma)
 %compute the log-density of multivariate normal
 %mu must be row vector
 [N,dim]=size(x);
 
 x=x-repmat(mu,N,1); %substract the mean

 inv_sigma=inv(sigma);%inverse covariance matrix
 
 
 left_x=repmat(x(:,1),1,dim);
 for i=2:dim
  left_x=horzcat(left_x,repmat(x(:,i),1,dim));
 end
 
 right_x=repmat(x,1,dim);
 f=(left_x.*right_x)*reshape(inv_sigma,[],1);%vec of (x-u)'inv(sigma)(x-u)
 f=f+ log(det(2*pi*sigma)); %normalizing factor
 f=-0.5*f;
 
