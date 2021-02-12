tstart=tic;

%parallel_pool=gcp('nocreate');

%% load parameters for DA
load 'state_par20210104_new.mat' % load model states and parameters
load 'state_par20210111_new.mat'
% load('beta_ratio.mat','beta_ratio') ; % load beta scenario
%load('Vdoses.mat','Vdoses');  % load vaccination
% load('timesV180.mat','timesV');  % load vaccination

%% select initial beta for each realization
beta_red_in=squeeze(beta_r(:,2,:)); % matrix with ensemble of initial betas (100 X 107)
time_red_in=time_b(2);

NSample=size(x0_real_out,1);
%NSample=10;
V.NSample=NSample;
x0_real=squeeze(x0_real_out(1:V.NSample,:,1)); % matrix with ensemble of initial states 
x0_real(:,end+1:end+2*V.n)=zeros(V.NSample,2*V.n); % add cumulative exposed and vaccinated

clear x0_real_out
%% set final time of the simulation 
time_final=timesV(end);
V.times=time_red_in:timesV(end)%+1;

V.gammaV = 0; %1 / (9 * 30); % loss of immunity from faccination
n=V.n;

%% loop on V.times
t0=tic;
disp('start loop on times and MC realizations')
disp(['current time: ',datestr(V.times(1)),...
    ' ; final time:', datestr(V.times(end))])
%% parameters to be passed in the parallel loop
ens_exposed=zeros(V.NSample,V.n,length(V.times)-1);
parfor cont_sample=1:V.NSample
    %disp('start real')
    % select current parameters and state variables for the realization
    Vt=V;
    PAR=PAR_real(cont_sample,:);
    
    Vt.x0=squeeze(x0_real(cont_sample,:))';
    Vt.ProvBetaInc=1;
    % define time_model if tm==1
    Vt.beta3_red_reg=(squeeze(beta_red_in(cont_sample,:))'*beta_ratio); %sequence of beta for current realization
    Vt.tbeta=Vt.times(1):(Vt.times(1)+length(beta_ratio));
    out_exposed=zeros(n,length(Vt.times)-1);
    for tm=2:length(Vt.times)
        Vt.time_model=Vt.times(tm-1):Vt.times(tm);
        
        % compute vaccination rate
        
        S=Vt.x0(1:n);
        E=Vt.x0(n+1:2*n);
        P=Vt.x0(2*n+1:3*n);
        A=Vt.x0(4*n+1:5*n);
        R=Vt.x0(7*n+1:8*n);
        
        ii=find(timesV==Vt.times(tm-1));
        
        v_rate=Vdoses(ii,:)'./(S+E+P+A+R);
        Vt.x0(1:n)=S-v_rate.*S;
        Vt.x0(12*n+1:13*n)=Vt.x0(12*n+1:13*n)+v_rate.*S;
        %% run model
        [x]=sepia_spatial_vac(PAR,Vt);
        
        %
        Vt.x0=x(2,:)';
        out_exposed(:,tm-1)=Vt.x0(11*n+1:12*n);
    end
    ens_exposed(cont_sample,:,:)=out_exposed;

    %disp('finish real')
end

ens_exposed_preprocess = sum(sum(ens_exposed,2),3);