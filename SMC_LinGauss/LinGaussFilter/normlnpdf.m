function f=normlnpdf(x,mu,sigma2)

f=-(x-mu).^2./(2*sigma2)-log(sigma2)/2-log(2*pi)/2;