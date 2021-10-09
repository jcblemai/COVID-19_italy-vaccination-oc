function LogLKNB=LogLKNB(x,m,omega)
%x data (MUST BE INTEGER)
%m average of the BN (vector same size of y)
%omega to inflate variance at omega*m

%p and r parameter of  the negative binomial distribution
x(x<1)=1;
omega1=omega+0*m;
%omega1(m<5)=omega1(m<5)./m(m<5)*5;
%omega1(abs(omega1-1)<0.01)=1.01;
r=zeros(size(m));
p=r;
r=m./(omega1-1);
p=1./omega1;
% increase variance for small model outputs
p(m<5)=m(m<5)./(5*omega1(m<5));
r(m<5)=m(m<5).*p(m<5)./(1-p(m<5));


LogLKNB = sum((gammaln(r+x) - gammaln(x+1) - gammaln(r) + r.*log(p) + x.*log1p(-p)));