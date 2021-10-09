function [LogPosterior,x,R0,R_t,e_t,beta_red_out]=sepia_spatial_DA(PAR,Vt,data)

%other parameters
n=Vt.n;  %number of nodes

% parameter values
deltaE=1/PAR(2);
deltaP=1/PAR(3);
sigma=PAR(4);
eta=1/PAR(5);
gammaH=1/PAR(6);
gammaI=1/PAR(6);
alphaI=0;
alphaH=1/40;%PAR(7);
gammaH=1/PAR(6);
epsilonP=Vt.epsilonP;
epsilonA=PAR(8);
rS=PAR(9);
%Deltat0=PAR(10);
epsilonI=PAR(11)*epsilonA;
omega=PAR(12);  %parameter r of NB distribution

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
on=ones(Vt.n,1);
beta3_reg=PAR(Vt.nPAR_model+1:Vt.nPAR_model+Vt.n_reg)';
beta3=beta3_reg(Vt.prov_IDreg).*Vt.ProvBetaInc;

if Vt.DA_par_flag
%     if Vt.time_model(1)>=Vt.tbeta(end)
%          beta3_old=Vt.beta3_reg_old;
%          %tbeta=[Vt.time_model(1),floor((Vt.time_model(1)+Vt.time_model(end))/2),Vt.time_model(end)];
%          tbeta=[Vt.time_model(1),floor((Vt.time_model(1)+Vt.time_model(end))/2),Vt.time_model(end)];
%          beta_p=[beta3_old,beta3, beta3];
%     else
%         if Vt.time_model(end)<=Vt.tbeta(end)
%             tbeta=Vt.tbeta;
%             beta_p=eval(Vt.beta_string);
%         else
%         %tbeta=[Vt.tbeta,floor((Vt.tbeta(end)+Vt.time_model(end))/2),Vt.time_model(end)];
%         %beta_p=[eval(Vt.beta_string),beta3,beta3];
%         tbeta=[Vt.tbeta,Vt.time_model(end)];
%         beta_p=[eval(Vt.beta_string),beta3];
%             
%         end
%     end
% else
%     tbeta=Vt.tbeta;
%     beta_p=eval(Vt.beta_string);
tbeta=[Vt.time_model(end)];
beta_p=[ beta3];
end
%% select mobility reduction from google (provices)
time_model=Vt.time_model;
for it=1:length(time_model)
    if ismember(time_model(it),Vt.mob_red_d)
        indt=find(Vt.mob_red_d==Vt.time_model(it));
        mob_red(:,it)=(1+Vt.mob_red_w(indt,:)'/100);
    elseif  time_model(it)<Vt.mob_red_d(1)
        mob_red(:,it)=(1+Vt.mob_red_w(1,:)'/100);
    elseif time_model(it)>Vt.mob_red_d(end)
        mob_red(:,it)=(1+Vt.mob_red_w(end,:)'/100);
    end
end
mob_red(mob_red<0.1)=0.1;
%calculate beta reduction matrix for output
%beta_p=beta_p.*mob_red;
%beta_red_out=interp1(tbeta,beta_p',time_model);
beta_red_out=repmat(beta_p,1,length(time_model));
%mob_red_out=interp1(tmob,mob',Vt.time_model);
beta_red_out=beta_red_out'.*(mob_red');


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
cumH=x(:,10*n+1:11*n)'; %cumulative cases

e_t=zeros(1,length(Vt.time_model));
R_t=zeros(1,length(Vt.time_model));

if Vt.time_model(1)<datenum(2020,02,24)
S=x(:,1:n)';
E=x(:,n+1:2*n)';
P=x(:,2*n+1:3*n)';
I=x(:,3*n+1:4*n)';
A=x(:,4*n+1:5*n)';
%Q=x(:,5*n+1:6*n)';
%H=x(:,6*n+1:7*n)';
R=x(:,7*n+1:8*n)';
%RR=x(:,8*n+1:9*n)'; %recorder recovered
%D=x(:,9*n+1:10*n)';


%compute endemicity
%tic

epsilon=zeros(n,1);
chiE=zeros(n,1);
chiP=zeros(n,1);
chiI=zeros(n,1);
chiA=zeros(n,1);
mu=(1/(82*365));
MI=diag(ones(n,1)); %mobility infected is null
%mob_red_out=interp1(tmob,mob',Vt.time_model);
mob_red_out=mob_red';
for i=1:length(time_model)
    betaP=diag(beta0*beta_red_out(i,:)*epsilonP);
    betaI=diag(beta0*beta_red_out(i,:)*epsilonI);
    betaA=diag(beta0*beta_red_out(i,:)*epsilonA);
%% 

   CS=(rS*Vt.p.*mob_red_out(i,:)).*Vt.q;       %compute only CS
    CS(1:n+1:end)=1-(sum(CS,2)-diag(CS));
    
    
    [R_t(i),e_t(i)]=SEPIAR_eigen_t(betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
        mu,deltaE,deltaP,sigma,gammaI,eta,alphaI,gammaA,...
        S(:,i),E(:,i),P(:,i),I(:,i),A(:,i),R(:,i),CS,CS,CS,MI,CS,CS,n,Vt.W,Vt.Wpi);
end
%R0=beta0*epsilonP/deltaP + epsilonI*sigma*beta0/(eta+gammaI+alphaI)+(1-sigma)*beta0*epsilonA/gammaA;
R0=R_t(1);
else
   R0=0; 
end

% compute posterior

%{'prov_Hnew';'reg_Hnew';'prov_Rnew';'reg_Rnew';'prov_Dnew';'reg_Dnew'};
LogPosterior.prov=0;
LogPosterior.reg=0;
LogPosterior.loc_prov=zeros(Vt.n_reg,1);
LogPosterior.loc_reg=zeros(Vt.n_reg,1);

if isempty(data.Date) 
    dd=0;          % for the projections, there aren't dates for computing LogLike
else
    dd=data.Date(end);%NOTE Vt.Date is updated in every DA window
end

if Vt.time_model(1)<dd
    if  Vt.LK(1)
        prov_cumH_sim=interp1(Vt.time_model,cumH',[data.Date(1)-1,data.Date])';
        prov_Hnew_sim=diff(prov_cumH_sim,1,2);
        
        LogPosterior.prov=LogPosterior.prov+LogLKNB(data.prov_Hnew(:),prov_Hnew_sim(:),omega);
        
        reg_Hnew_sim=Vt.prov2reg*prov_Hnew_sim;
        LogPosterior.reg =LogPosterior.reg+LogLKNB(data.reg_Hnew(:),reg_Hnew_sim(:),omega);
        
        for reg=1:Vt.n_reg
            prov_data=data.prov_Hnew(Vt.prov_IDreg==reg,:);
            prov_sim=prov_Hnew_sim(Vt.prov_IDreg==reg,:);
            
            LogPosterior.loc_prov(reg)=LogPosterior.loc_prov(reg)+LogLKNB(prov_data(:),prov_sim(:),omega);
            
            LogPosterior.loc_reg(reg)=LogPosterior.loc_reg(reg)+LogLKNB(data.reg_Hnew(reg,:),reg_Hnew_sim(reg,:),omega);
        end
        
    else
        error('updated computation of loglike for other measures') 
    end
    %%%% OTHER OPTIONS HAVE NOT BEEN UPDATED FOR data assimilation
%     if Vt.LK(2)
%         reg_cumH_sim=Vt.prov2reg*(interp1(Vt.time_model,cumH',Vt.Date)');
%         reg_Hnew_sim=diff(reg_cumH_sim,1,2);
%         LogPosterior=LogPosterior+LogLKNB(Vt.reg_Hnew(:),reg_Hnew_sim(:),omega);
%     end
%     if Vt.LK(3)
%         prov_RR_sim=interp1(Vt.time_model,RR',Vt.Date)';
%         prov_RRnew_sim=diff(prov_RR_sim,1,2);
%         LogPosterior=LogPosterior+LogLKNB(Vt.prov_Rnew(:),prov_RRnew_sim(:),omega);
%     end
%     if Vt.LK(4)
%         reg_RR_sim=Vt.prov2reg*(interp1(Vt.time_model,RR',Vt.Date)');
%         reg_RRnew_sim=diff(reg_RR_sim,1,2);
%         LogPosterior=LogPosterior+logLKNB(Vt.reg_Rnew(:),reg_RRnew_sim(:),omega);
%     end
%     if Vt.LK(5)
%         prov_D_sim=interp1(Vt.time_model,D',Vt.Date)';
%         prov_Dnew_sim=diff(prov_D_sim,1,2);
%         LogPosterior=LogPosterior+LogLKNB(Vt.prov_Dnew(:),prov_Dnew_sim(:),omega);
%     end
%     if Vt.LK(6)
%         reg_D_sim=Vt.prov2reg*(interp1(Vt.time_model,D',Vt.Date)');
%         reg_Dnew_sim=diff(reg_D_sim,1,2);
%         LogPosterior=LogPosterior+LogLKNB(Vt.reg_Dnew(:),reg_Dnew_sim(:),omega);
%     end
%    if sum(Vt.prior)>0
%        LogPosterior=LogPosterior-sum((PAR(Vt.prior)'-Vt.prior_par(Vt.prior,1)).^2./(2*Vt.prior_par(Vt.prior,2).^2));
%    end
end

%%R0=beta0*epsilonP/deltaP + epsilonI*sigma*beta0/(eta+gammaI+alphaI)+(1-sigma)*beta0*epsilonA/gammaA;

%if Vt.show_R0
%    R0
%    R0P=beta0*epsilonP/deltaP
%    R0I=epsilonI*sigma*beta0/(eta+gammaI+alphaI)
%    R0A=(1-sigma)*beta0*epsilonA/gammaA
%end

%**********************************************************************************
%NESTED FUNCTION
%**********************************************************************************
    function dxdt=eqs(t,x)
        
        %compute beta reduction
        %if t<=tbeta(1)
        %    beta_red=beta_p(:,1);
        %elseif t>=tbeta(end)
        %    beta_red=beta_p(:,end);
        %else
        %    ii=find(tbeta<t,1,'last');
        %    m=(beta_p(:,ii+1)-beta_p(:,ii))/(tbeta(ii+1)-tbeta(ii));
        %    beta_red=(t-tbeta(ii))*m+beta_p(:,ii);
        %end
        indm=find(time_model>=t,1,'first');
        beta_red=beta_p.*mob_red(:,indm);
        %compute mobility reduction
%         if t<tmob(1)
%             mob_red=mob(:,1);
%         elseif t>tmob(end)
%             mob_red=mob(:,end);
%         else
%             ii=find(tmob<t,1,'first');
%             m=(mob(:,ii+1)-mob(:,ii))/(tmob(ii+1)-tmob(ii));
%             mob_red=(t-tmob(ii))*m+mob(:,ii);
%         end
        
        CS=(rS*Vp.*mob_red(:,indm)).*Vq;       %compute only CS
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
        
        %Nmob=CS'*S+CE'*E+CP'*P+CI'*I+CA'*A+CR'*R;
        % FoI=CS*((betaP*CP'*P+betaI*CI'*I+betaA*CA'*A)./Nmob).*S;
        
        %simplified version with only one mobility matrix for S E P A R and
        %I do not move
        Nmob=CS'*(S+E+P+R+A)+I;
        %assumption. class I do not move
        FoI=CS*((CS'*(beta0*beta_red.*(epsilonP*P+epsilonA*A))+epsilonI*beta0*beta_red.*I)./Nmob).*S;
        
        dSdt=-FoI;
        dEdt=FoI-deltaE*E;
        dPdt=deltaE*E-deltaP*P;
        dIdt=sigma*deltaP*P-(eta+gammaI+alphaI)*I;
        dAdt=(1-sigma)*deltaP*P-gammaA*A;
        dQdt=(1-VfracH)*eta*I-gammaQ.*Q;
        dHdt=VfracH*eta*I-(gammaH+alphaH)*H;
        dRdt=gammaI*I+gammaA*A+gammaH*H+gammaQ*Q;
        dRRdt=gammaH*H;
        dDdt=alphaH*H;
        dcumHdt=VfracH*eta*I;
        
        dxdt=[dSdt; dEdt; dPdt; dIdt; dAdt; dQdt; dHdt; dRdt; dRRdt; dDdt; dcumHdt];
        return
    end
return
end
