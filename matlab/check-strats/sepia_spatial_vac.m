function [x_out]=sepia_spatial_mod(PAR,Vt)

%other parameters
n=Vt.n;  %number of nodes

% parameter values
deltaE=1/PAR(2);
deltaP=1/PAR(3);
sigma=PAR(4);
eta=1/PAR(5);
gammaI=1/PAR(6);
alphaI=1/PAR(7);
alphaH=1/PAR(7);
gammaH=1/PAR(6);
epsilonP=Vt.epsilonP;
epsilonA=PAR(8);
rS=PAR(9);
%Deltat0=PAR(10);
epsilonI=PAR(11)*epsilonA;
omega=PAR(12);  %parameter r of NB distribution
gammaV=Vt.gammaV;

%beta1=PAR(13);
%beta2=PAR(14);
%mean_beta3=PAR(15); % no need for these parameters here 
%std_beta3=PAR(16);  % no need for these parameters here


gammaQ=Vt.gammaQgammaH*gammaH;
gammaA=Vt.gammaAgammaQ*gammaQ;


beta0=PAR(1)/(1/deltaP + epsilonI*sigma/(gammaI + alphaI + eta) + epsilonA*(1-sigma)/gammaA);

%variables for mobility reduction matrix
mob=Vt.mob;
tmob=Vt.tmob;

% variables for beta reduction

tbeta=Vt.tbeta;
beta_p=Vt.beta3_red_reg;
%mob_red_out=interp1(tmob,mob',Vt.time_model);

Vp=Vt.p;
Vq=Vt.q;

VfracH=Vt.fracH;

% simulation
%options=odeset('NonNegative',Vt.neqs);
%options=odeset('RelTol',1e-7,'AbsTol',1e-8);
options=odeset;

[t,x]=ode45(@eqs,Vt.time_model,Vt.x0,options);
%sol=ode45(@eqs,Vt.time_model,Vt.x0,options);

%postprocessing
if length(Vt.time_model)==2
il=[1,size(x,1)];
else
    il=1:length(Vt.time_model);
end

x_out=x(il,:);

%**********************************************************************************
%NESTED FUNCTION
%**********************************************************************************
    function dxdt=eqs(t,x)
        
        %compute beta reduction
        if t<=tbeta(1)
            beta_red=beta_p(:,1);
        elseif t>=tbeta(end)
            beta_red=beta_p(:,end);
        else
            ii=find(tbeta<t,1,'last');
            m=(beta_p(:,ii+1)-beta_p(:,ii))/(tbeta(ii+1)-tbeta(ii));
            beta_red=(t-tbeta(ii))*m+beta_p(:,ii);
        end
        
        %compute mobility reduction
        if t<tmob(1)
            mob_red=mob(:,1);
        elseif t>tmob(end)
            mob_red=mob(:,end);
        else
            ii=find(tmob<t,1,'first');
            m=(mob(:,ii+1)-mob(:,ii))/(tmob(ii+1)-tmob(ii));
            mob_red=(t-tmob(ii))*m+mob(:,ii);
        end
        
        CS=(rS*Vp.*mob_red).*Vq;       %compute only CS
        %       CE=Vt.p.*Vt.q*rE; CP=Vt.p.*Vt.q*rP;
        %       CI=Vt.p.*Vt.q*rI; CA=Vt.p.*Vt.q*rA; CR=Vt.p.*Vt.q*rR;
        %
        %set diagonal imposing row-stocasticity
        CS(1:n+1:end)=1-(sum(CS,2)-diag(CS));
        %CS(1:n+1:end)=(1-Vt.p')+(1-rS)*Vt.p'+rS*Vt.p'.*Vt.q(1:n+1:end);
        %         CE(1:n+1:end)=(1-Vt.p')+(1-rE)*Vt.p'+rE*Vt.p'.*Vt.q(1:n+1:end);
        %         CP(1:n+1:end)=(1-Vt.p')+(1-rP)*Vt.p'+rP*Vt.p'.*Vt.q(1:n+1:end);
        %         CI(1:n+1:end)=(1-Vt.p')+(1-rI)*Vt.p'+rI*Vt.p'.*Vt.q(1:n+1:end);
        %         CA(1:n+1:end)=(1-Vt.p')+(1-rA)*Vt.p'+rA*Vt.p'.*Vt.q(1:n+1:end);
        %         CR(1:n+1:end)=(1-Vt.p')+(1-rR)*Vt.p'+rR*Vt.p'.*Vt.q(1:n+1:end);
        %
        
        S=x(1:n);
        E=x(n+1:2*n);
        P=x(2*n+1:3 *n);
        I=x(3*n+1:4*n);
        A=x(4*n+1:5*n);
        Q=x(5*n+1:6*n);
        H=x(6*n+1:7*n);
        R=x(7*n+1:8*n);
%         RR=x(8*n+1:9*n); %recorder recovered
%         D=x(9*n+1:10*n);
%         cumH=x(10*n+1:11*n); %cumulative cases
%   
%        cumE=x(11*n+1:12*n)
        VV=x(12*n+1:13*n);
        
        %Nmob=CS'*S+CE'*E+CP'*P+CI'*I+CA'*A+CR'*R;
        % FoI=CS*((betaP*CP'*P+betaI*CI'*I+betaA*CA'*A)./Nmob).*S;
        
        %simplified version with only one mobility matrix for S E P A R and
        %I do not move
        Nmob=CS'*(S+E+P+R+A)+I;
        %assumption. class I do not move
        FoI=CS*((CS'*(beta0*beta_red.*(epsilonP*P+epsilonA*A))+epsilonI*beta0*beta_red.*I)./Nmob).*S;
        
        dSdt=-FoI;%+gammaV*VV;
        dEdt=FoI-deltaE*E;
        dPdt=deltaE*E-deltaP*P;
        dIdt=sigma*deltaP*P-(eta+gammaI)*I; % +alphaI
        dAdt=(1-sigma)*deltaP*P-gammaA*A;
        dQdt=(1-VfracH)*eta*I-gammaQ.*Q;
        dHdt=VfracH*eta*I-(gammaH+alphaH)*H;
        dRdt=gammaI*I+gammaA*A+gammaH*H+gammaQ*Q;
        dRRdt=gammaH*H;
        dDdt=alphaH*H;
        dcumHdt=VfracH*eta*I;
        dVdt=-gammaV*VV;
        dcumE=FoI;
        
        dxdt=[dSdt; dEdt; dPdt; dIdt; dAdt; dQdt; dHdt; dRdt; dRRdt; dDdt; dcumHdt;dcumE;dVdt;];
        return
    end
return
end
