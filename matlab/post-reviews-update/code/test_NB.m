x=round(0:400)



m=2
omega=5
r=m/(omega-1);
p=1/omega;

figure
plot(x,nbinpdf(x,r,p))


return
%%
N_real=1e6

m=20
omega=5
r=m/(omega-1);
p=1/omega;
ave=((1-p)*r)/(p)
variance=((1-p)*r)/((p)^2)
X=nbinrnd(r,p,N_real,1);
CI_X=quantile(X,[0.025 0.5 0.975]);

m=30;
omega=5;
r=m/(omega-1);
p=1/omega;
ave=((1-p)*r)/(p)
variance=((1-p)*r)/((p)^2)
Y=nbinrnd(r,p,N_real,1);
CI_Y=quantile(Y,[0.025 0.5  0.975]);


Z=X+Y;

CI_Z=quantile(Z,[0.025 0.5  0.975])
CI_X+CI_Z


%%
ave=mean(X+Y)
variance=var(X+Y)
