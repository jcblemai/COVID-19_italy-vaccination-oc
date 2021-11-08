close all
tstart=tic;

%parallel_pool=gcp('nocreate');

%% load parameters for DA
%load 'state_par20210104_new.mat' % load model states and parameters
%load 'state_par20210111_new.mat'
%load 'input_20211012/google_data_and_ages.mat'
%load 'input_20211101/state_par2020-04-30.mat'
% load 'input_20211101/state_par2020-12-31.mat'
load 'input_20211101/state_par2021-01-04.mat'

%load('beta_ratio.mat','beta_ratio') ; % load beta scenario
%load('Vdoses.mat','Vdoses');  % load vaccination
%load('timesV180.mat','timesV');  % load vaccination

%beta_ratio=beta_ratio(end:-1:1)/1.5;

% NOTE: timesV and Vdoses have different numbers of dates

flag_age_classes=1;% 1 to compute also age classes; 0 otherwise

%% temp
load input_20211101/beta_red_out2021-01-30 % load all values of beta

NSample=size(x0_real,1);
%NSample=10;
V.NSample=NSample;
x0_real=squeeze(x0_real(1:V.NSample,:,1)); % matrix with ensemble of initial states
x0_real(:,end+1:end+2*V.n)=zeros(V.NSample,2*V.n); % add cumulative exposed and vaccinated

%% set final time of the simulation
time_final=timesV(end);
%time_final=datenum(2020,06,01);

V.times=time_x0:(time_final+1);%+1;

%% select beta values 

%indb=find(time_b>=timesV(1),1,'first');
%beta_red_in=squeeze(beta_r(:,indb,:));  % matrix with ensemble of initial betas (100 X 107)

indb=find(time_b<=timesV(1),1,'last');
beta_red_in=squeeze(beta_r(:,indb,:));  % matrix with ensemble of initial betas (100 X 107)


%% google mobility data
%indm=find(V.mob_red_d>=timesV(1),1,'first');
%mob_red=(1+V.mob_red_w(indm+(1:length(beta_ratio)),:)/100);
%beta_ratio=(diag(beta_ratio)*mob_red)';

clear x0_real_out


V.gammaV = 0; %1 / (9 * 30); % loss of immunity from faccination
n=V.n;

%% loop on V.times
t0=tic;
disp('start loop on times and MC realizations')
disp(['current time: ',datestr(V.times(1)),...
    ' ; final time:', datestr(V.times(end))])
%% parameters to be passed in the parallel loop
ens_exposed=zeros(V.NSample,V.n,length(V.times)-1);
if flag_age_classes
    ens_agec_E=zeros(V.NSample,V.n,V.nc,length(V.times)-1);
    ens_agec_S=zeros(V.NSample,V.n,V.nc,length(V.times)-1);
    ens_agec_D=zeros(V.NSample,V.n,V.nc,length(V.times)-1);
end

%%
nc=V.nc;
timesm=V.times;
mob_red_d=V.mob_red_d;
mob_red_w=V.mob_red_w;
gammaQgammaH=V.gammaQgammaH;
gammaAgammaQ=V.gammaAgammaQ;
Vp=V.p;
Vq=V.q;
VfracH=V.fracH;
epsilonP=V.epsilonP;
gammaV=V.gammaV;
agec_coef_pos=V.agec_coef_pos;
agec_letality=V.agec_letality;
n=V.n;
%%
parfor cont_sample=1:V.NSample
    
    Vt=[];
    %disp('start real')
    % select current parameters and state variables for the realization
    %Vt=V;
    PAR=PAR_real(cont_sample,:);
    
    x0=squeeze(x0_real(cont_sample,:))';
    %
    %Vt.beta3_red_reg=repmat(beta_red_in(cont_sample,:)',1,size(beta_ratio,2)).*beta_ratio; %sequence of beta for current realization
   
    %Vt.beta3_red_reg=squeeze(beta3_out(cont_sample,:,:));
    
    
    %Vt.tbeta=Vt.times(1):(Vt.times(1)+length(beta_ratio));
    if flag_age_classes
        agec_S1=squeeze(agec_S1_real(cont_sample,:,:));
        agec_E=zeros(n,nc,length(timesm)-1);
        agec_D=zeros(n,nc,length(timesm)-1);
    else
        
    end
    out_exposed=zeros(n,length(timesm)-1);
    for tm=2:length(timesm)
        Vt.time_model=timesm(tm-1):timesm(tm);
        
        indm=find(mob_red_d>=timesm(tm-1),1,'first');
        Vt.mob_red=squeeze((1+mob_red_w(indm,:)/100));

%      
        
        indt=find(timesV==timesm(tm-1));
        
        if length(indt)>0
           % disp('vaccination')
            % compute vaccination rate
            
            S=x0(1:n);
            E=x0(n+1:2*n);
            P=x0(2*n+1:3*n);
            A=x0(4*n+1:5*n);
            R=x0(7*n+1:8*n);
            if size(Vdoses)==3
                Vdoses1=squeeze(sum(squeeze(Vdoses(indt,:,:),2))); % sum on age classes
                Vdoses_agec=squeeze(Vdoses(indt,:,:));
            else
                Vdoses1=Vdoses(indt,:);
                if flag_age_classes
                    Vdoses_agec=diag(Vdoses1)* (agec_S1./sum(agec_S1,2)); % distribute vaccines accordingly to susceptibles
                end 
            end
          
            % vaccinate individuals per node
            v_rate=Vdoses1'./(S+E+P+A+R);
            
            
            x0(1:n)=S-v_rate.*S;
            x0(n+1:2*n)=E-v_rate.*E;
            x0(2*n+1:3*n)=P-v_rate.*P;
            x0(4*n+1:5*n)=A-v_rate.*A;
            x0(7*n+1:8*n)=R-v_rate.*R;
            
            x0(x0<0)=0;
            
            x0(12*n+1:13*n)=x0(12*n+1:13*n)+Vdoses1';
            
            if flag_age_classes
                
                % vaccinate individuals per age class
                v_rateS=v_rate.*S;
                Vdoses_agec_S=diag(v_rateS)*(Vdoses_agec./sum(Vdoses_agec,2));
                agec_S1=agec_S1-Vdoses_agec_S;
                agec_S1(agec_S1<0)=0;
            end
            
            %update beta3
           
           Vt.beta3_red=beta_ratio(indt)*squeeze(beta_red_in(cont_sample,:))'.*Vt.mob_red' ;   
        else
            Vdoses1=zeros(n,1);
            Vdoses_agec=zeros(n,1);

            indb=find(time_red_out>=timesm(tm-1),1,'first');
            Vt.beta3_red=squeeze(beta3_out(cont_sample,indb,:)).*Vt.mob_red';   
        end
        Vt.x0=x0;
        %Vt.tbeta=Vt.time_model(1);
        Vt.gammaQgammaH=gammaQgammaH;
        Vt.gammaAgammaQ=gammaAgammaQ;
        Vt.p=Vp;
        Vt.q=Vq;
        Vt.fracH=VfracH;
        Vt.epsilonP=epsilonP;
        Vt.gammaV=gammaV;
        Vt.n=n;
        
              
        %% run model
        [x]=sepia_spatial_vac(PAR,Vt);
        
        x0=x(2,:)';  %ic for next step
        
        % select exposed
        cumE=x0(11*n+1:12*n);
        out_exposed(:,tm-1)=cumE;
        dE=x(2,11*n+1:12*n)' -x(1,11*n+1:12*n)'; 
        
        
        
        %% compute cases per age class per province
        if flag_age_classes
            
            % update probability of becoming infected based on the susceptibles
            prov_S=x(1,1:n);
            
            prob_pos_prov=diag(1./sum(agec_S1,2))*(agec_coef_pos.*agec_S1);
            
            ss1=sum(prob_pos_prov,2);
            prob_pos_prov=prob_pos_prov./ss1; % riscale probabilities so that sum to 1
            agec_E(:,:,tm-1)= diag(dE)*prob_pos_prov; % distribute new cases per age class
            
            agec_S1=agec_S1-agec_E(:,:,tm-1);
            agec_S1(agec_S1<0)=0;
            %agec_S(:,:,tm-1)=agec_S1;
            dead=repmat(agec_letality,n,1).*squeeze(agec_E(:,:,tm-1));
            agec_D(:,:,tm-1)=dead;
        end%
    end
    ens_exposed(cont_sample,:,:)=out_exposed;
    if flag_age_classes
        
        ens_agec_E(cont_sample,:,:,:)=agec_E;
        %ens_agec_S(cont_sample,:,:,:)=agec_S;
        ens_agec_D(cont_sample,:,:,:)=agec_D; % note that times to death should be traslated of 14days
    end%disp('finish real')
end

ens_exposed_preprocess = sum(sum(ens_exposed,2),3);

return

figure
plot(timesm(3:end),diff(squeeze(sum(ens_exposed,2))')' )
datetick('x','dd-mm')
title('Ensemble exposed')

if flag_age_classes
EE=squeeze(sum(ens_agec_E,2));
EEq=squeeze(quantile(EE,0.5,1));
cEEq=cumsum(EEq);
figure
plot(timesm(2:end),cEEq);
datetick('x','dd-mm')
legend('0-19','20-39','40-59','60-79','80+')
title('Exposed (cumulative per age class, median)')

DD=squeeze(sum(ens_agec_D,2));
DDq=squeeze(quantile(DD,0.5,1));
cDDq=cumsum(DDq);
figure
plot(timesm(2:end)+V.agec_time_death,cDDq);
legend('0-19','20-39','40-59','60-79','80+')
datetick('x','dd-mm')
title('Deaths (cumulative per age class, median)')
end
