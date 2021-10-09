
function [L_naz1,L_reg1,L_naz_for1,L_reg_for1]=function_main(NSample,NIter,NLagDays,NObs,V,Posterior_sample)

V.NSample=NSample;
V.NIter=NIter;
V.NLagDays=NLagDays;
V.NObs=NObs;

V.filename=['NS',num2str(NSample),...
    '_NI',num2str(NIter),'_NL',num2str(NLagDays),'_NO',num2str(NObs)];
%if ~exist(V.folder_name, 'dir')
%       mkdir(V.folder_name)
%end

clear NSample NIter NLagDays NObs


%% initialize V.times for DA
if V.DA_flag
    if V.time_first_DA>V.time_final
        error('First assimilation after the end of simulation, check dates for DA')
    end
    
    if V.time_final>V.data.Date(end)-V.NLagDays
        V.times=V.time_first_DA:1:(V.data.Date(end)-V.NLagDays-V.NObs); % forecast will be up to time_final
    else
        V.times=V.time_first_DA:1:(V.time_final-V.NLagDays-V.NObs);
    end
else
    V.times=V.time_final; % no assimilation
    V.forecast_days=0;
    V.NSample=1;
end

%% Create Posteriror Sample and select MC realizations
PAR_real=zeros(V.NSample,size(Posterior_sample,2));

%mean_beta3=PAR(15);
%std_beta3=PAR(16);

vec1=zeros(V.NSample,V.n_reg);
rs=rand(V.NSample,1);
DA_par_flag=V.DA_par_flag;
min_truncate=V.min_truncate;
max_truncate=V.max_truncate;
std_beta3i=V.std_beta3i;
n_reg=V.n_reg;
epsilonP=V.epsilonP;
parfor cont_sample=1:V.NSample
%for cont_sample=1:V.NSample
 
    ind=ceil(rs(cont_sample)*size(Posterior_sample,1));
    PAR_real(cont_sample,:)=Posterior_sample(ind,:);
   
    if DA_par_flag
        % initialize beta3 per region
        %beta1=ps(13);
        %beta2=ps(14);
        gt=makedist('Normal','mu',1,'sigma',std_beta3i);
        gt=truncate(gt,min_truncate,max_truncate);
        %gt=makedist('Normal','mu',1,'sigma',V.std_beta3i);
        %gt=truncate(gt,V.min_truncate,V.max_truncate);
        beta3=random(gt,1,n_reg);
        %vec1(cont_sample,:)=beta1*beta2*beta3;
        vec1(cont_sample,:)=beta3;
    end
end
PAR_real(:,(V.nPAR_model+1):(V.nPAR_model+V.n_reg))=vec1;
%DA.PAR_real=PAR_real;
clear i std_beta3 Posterior_sample vec1 i rs

%% initialize log weights for DA
DA.logw_real=-ones(V.NSample,1)*log(V.NSample);

DA.x_real = [];
DA.t_real = [];
DA.beta_red_real = [];
DA.LogL_naz_reg_real = [];
DA.LogL_loc_reg_real = [];
DA.Rt_real = [];
DA.Et_real = [];
DA.countryR_q=[];
DA.countryE_q=[];
beta_red_out=[];
t_red_out=[];
%L_reg1=[];
%L_naz1=[];
L_reg_for1=[];
L_naz_for1=[];


count=0;
count_for=0;

%% initial conditions
% define initial conditions with seeding
x0_real=repmat(V.x0',V.NSample,1);
x0_real(:,V.n+V.seeding)=10.^(PAR_real(:,V.nPAR_model+V.n_reg+1:V.nPAR_model+V.n_reg+length(V.seeding)));
time_seed_real=(V.data.Date(1)-(PAR_real(:,10)));

%% loop on V.times
t0=tic;
disp('start loop on times and MC realizations')
disp(['current time: ',datestr(V.data.Date(1)),...
    ' ; final time:', datestr(V.times(end))])
%% parameters to be passed in the parallel loop
n=V.n;
gammaQgammaH=V.gammaQgammaH;
gammaAgammaQ=V.gammaAgammaQ;
mob=V.mob;
tmob=V.tmob;
nPAR_model=V.nPAR_model;
n_reg=V.n_reg;
prov_IDreg=V.prov_IDreg;
DA_par_flag=V.DA_par_flag;
tbeta=V.tbeta;
beta_string=V.beta_string;
Vp=V.p;
Vq=V.q;
fracH=V.fracH;
prov2reg=V.prov2reg;
data_Date1=V.data.Date(1);
time_first_DA=V.time_first_DA;
VWpi=V.Wpi; % pseudoinverse of W
VW=V.W;
LK=V.LK;

mob_red_w=V.mob_red_w;
mob_red_d=V.mob_red_d;
agec=V.agec;

for tm=1:length(V.times)
    timestm=V.times(tm);

    %% select data for the current DA window
    % for the computation of the LogPosterior in the assimilation window
    if V.DA_flag
        if tm==1
            % indices for selecting the correct data
            i1=find(V.data.Date(2:end)<=V.times(1),1,'first')+1;
            i2=find(V.data.Date(2:end)<=V.times(1)+V.NLagDays+V.NObs,1,'last')+1;
            if length(i1)>0
            time_model=min(V.data.Date(i1),(V.times(1)-1)):(V.times(1)+V.NLagDays+V.NObs);
            else
                time_model=V.times(1)-1:V.times(1);
            end
            ind_out=length([0,time_model]);
            
            beta3_reg_old_real=ones(V.NSample,V.n);
        else
            time_model=V.times(tm-1):(V.times(tm-1)+V.NLagDays+V.NObs);%+V.forecast_days);            
            ind_out=length(time_model);
            if length(time_model)<=2
                error('check V.time_model')
            end
            i1=find(V.data.Date(2:end)>(time_model(end)-V.NObs),1,'first')+1;
            i2=find(V.data.Date(2:end)<=(time_model(end)),1,'last')+1;
            
            disp(' ')
            disp(['time:', datestr(V.times(tm-1)),...
                ' ; progress % ',num2str((V.times(tm-1)-V.times(1))/(V.times(end)-V.times(1)))])
            %' ; final t: ', datestr(V.times(end)), ...
            
  
            if tm==2
                beta3_reg_old_real=beta_red_real(:,end,:);
                %time_red_out=time_model(1);
                %time_x0=time_model(1);

                %beta_red_out=beta3_reg_old_real;
                %save(['results/state_par',datestr(time_model(1))],'V','PAR_real','time_red_out','beta_red_out','x0_real','time_x0')
            else
                beta3_reg_old_real=beta_red_real(:,2,:);
                %time_red_out(end+1)=time_model(1);
                %beta_red_out(:,end+1,:)=beta3_reg_old_real;
                
                %if ismember(time_model(1),[datenum(2020,11,1:5:30)])
                %     time_x0=time_model(1);

                %    save(['results/state_par',datestr(time_model(1))],'V','PAR_real','time_red_out','beta_red_out','x0_real','time_x0')
                %end
            end
            
            
        end
    else
        i1=2;
        i2=find(V.data.Date(2:end)<=V.times(1),1,'last')+1;
        ind_out=length([0,data.Date(1):(V.times(tm))]);
    end
    
    x_real=zeros(V.NSample,ind_out,V.neqs);
    t_real=zeros(V.NSample,ind_out);
    Rt_real=zeros(V.NSample,ind_out);
    Et_real=zeros(V.NSample,ind_out);
    beta_red_real=zeros(V.NSample,ind_out,V.n);
    R0_real=zeros(V.NSample,1);
    %LogL_naz_prov_real=zeros(V.NSample,1);
    LogL_naz_reg_real =zeros(V.NSample,1);
    %LogL_loc_prov_real=zeros(V.NSample,V.n_reg);
    LogL_loc_reg_real =zeros(V.NSample,V.n_reg);
    
    data.Date=V.data.Date(i1:i2);
    data_Date1=data.Date;
    
    if length(data.Date)>0
            data.prov_Hnew= V.data.prov_Hnew(:,(i1:i2)-1);
            data.reg_Hnew = V.data.reg_Hnew (:,(i1:i2)-1);
            data.reg_Rnew = V.data.reg_Rnew (:,(i1:i2)-1);
            data.reg_Dnew = V.data.reg_Dnew (:,(i1:i2)-1);
    end
    
    clear i1 i2 ind_out
    
    %% start loop on realizations
    t1=tic;
    
    
    NIter=1;
    if length(data_Date1)>V.NObs/2  % first model window: no DA
        if tm==1 
            NIter=4; 
        else
            NIter=V.NIter;
        end
    end
    
    for iterDA=1:NIter % iterations for parameter update
        %t2=tic;
        %disp(['DA iteration  ',num2str(iterDA)])
        parfor cont_sample=1:V.NSample
        %for cont_sample=1:V.NSample
            %for cont_sample=1:V.NSample
            Vt=struct;
            Vt.n=n;
            Vt.gammaQgammaH=gammaQgammaH;
            Vt.gammaAgammaQ=gammaAgammaQ;
            Vt.mob=mob;
            Vt.tmob=tmob;
            Vt.nPAR_model=nPAR_model;
            Vt.n_reg=n_reg;
            Vt.prov_IDreg=prov_IDreg;
            Vt.DA_par_flag=DA_par_flag;
            Vt.tbeta=tbeta;
            Vt.beta_string=beta_string;
            Vt.p=Vp;
            Vt.q=Vq;
            Vt.fracH=fracH;
            Vt.prov2reg=prov2reg;
            Vt.data_Date1=data_Date1;
            Vt.time_first_DA=time_first_DA;
            Vt.timestm=timestm;
            Vt.epsilonP=epsilonP;
            Vt.Wpi=VWpi; % pseudoinverse of W
            Vt.W=VW;
            Vt.LK=LK;
            Vt.mob_red_w=mob_red_w;
            Vt.mob_red_d=mob_red_d;
            % select current parameters and state variables for the realization
            PAR=PAR_real(cont_sample,:);
            Vt.x0=x0_real(cont_sample,:);
            Vt.ProvBetaInc=1;
            % define time_model if tm==1
            if tm==1
                % V.time_model is defined for each random realization
                % define time_model and add initial seeding in the exposed class
                if ~isempty(Vt.data_Date1)
                Vt.time_model=[time_seed_real(cont_sample),time_model];
                else
                    Vt.time_model=[time_seed_real(cont_sample),Vt.timestm-1,Vt.timestm];
                end
            else
                Vt.time_model=time_model;
            end
            
            if Vt.DA_par_flag
                Vt.beta3_reg_old=beta3_reg_old_real(cont_sample,:)';
            end
            
            %% run model
            [LogPosterior,x,R0,R_t,e_t,beta_red]=sepia_spatial_DA(PAR,Vt,data);
            
            %RR=x(1:i1,8*V.n+1:9*V.n)'; %recorder recovered
            %D=x(1:i1,9*V.n+1:10*V.n)';
            %prov_cumH=x(1:i1,10*V.n+1:11*V.n)'; %cumulative cases
            
            
            x_real(cont_sample,:,:)=x; 
            t_real(cont_sample,:)  =Vt.time_model;
            Rt_real(cont_sample,:) =R_t; 
            Et_real(cont_sample,:) =e_t;
            beta_red_real(cont_sample,:,:)=beta_red;
   
            %V.LogL_real.prov(cont_sample,count:(count+ind_out-1))=LogPosterior.prov*ones(1,ind_out);
            
            if tm==1
            R0_real(cont_sample)=R0;
            end
            %LogL_naz_prov_real(cont_sample)  =LogPosterior.prov;
            LogL_naz_reg_real (cont_sample)  =LogPosterior.reg;
            %LogL_loc_prov_real(cont_sample,:)=LogPosterior.loc_prov;
            LogL_loc_reg_real (cont_sample,:)=LogPosterior.loc_reg;
        end
        
        clear beta_red cont_sample ind_out x LogPosterior R0
        
        ind_out=size(x_real,2);
        
        %% extract variables to plot and compute quantiles
        if tm==1
            V.R0_real=R0_real;
        end
        
        DaTimes=mean(t_real,1);
        
        [regHnew_quantile,regR_quantile,regSfrac_quantile,countryHnew_quantile,countryR_quantile,countryE_quantile,regBeta_quantile]=extract_quantiles(PAR_real,V,x_real,DaTimes,beta_red_real,Rt_real,Et_real);
        %DA.x_real(:,(count+1):(count+ind_out),:)=x_real;
        DA.regHnew_q(:,(count+1):(count+ind_out-1),:)=regHnew_quantile;
        
        DA.countryHnew_q((count+1):(count+ind_out-1),:)=countryHnew_quantile;
        DA.times((count+1):(count+ind_out-1)) = DaTimes(2:end);
        
        %DA.beta_red_real(:,(count+1):(count+ind_out),:)=beta_red_real;
        DA.regR_q(:,(count+1):(count+ind_out-1),:)=regR_quantile(:,2:end,:);
        DA.regBeta_q(:,(count+1):(count+ind_out-1),:)=regBeta_quantile(:,2:end,:);
        
        
        DA.regSfrac_q(:,(count+1):(count+ind_out-1),:)=regSfrac_quantile(:,2:end,:);
        DA.LogL_naz_reg_real(:,(count+1))  =LogL_naz_reg_real;
        DA.LogL_loc_reg_real(:,(count+1),:)=LogL_loc_reg_real;
        
        DA.countryR_q((count+1):(count+ind_out-1),:)=countryR_quantile;
        DA.countryE_q((count+1):(count+ind_out-1),:)=countryE_quantile;
        
        
        time_red_out((count+1):(count+ind_out-1)) = DaTimes(2:end);
        beta_red_out(:,(count+1):(count+ind_out-1),:)=beta_red_real(:,2:end,:);
        
        
        %disp(['DA iteration  ',num2str(iterDA)])
        
        %disp(['DA iteration  ',num2str(iterDA) ', sim. time:',num2str(toc(t2)),' s'])
        
        %if   V.show_fig
        %    V.nfig_iter=V.nfig_iter+1;
        %    V.nfig=V.nfig+1;
        %    plot_DA_results_iter(DA,V,iterDA,time_model)
        %else
        %    if tm==length(V.times)&& iterDA==NIter
        %        plot_DA_results_iter(DA,V,iterDA,time_model)
        %    end
        %end
        
        %% update step of DA
        x0_real_temp=squeeze(x_real(:,1,:)); %state variables to be updated
        if V.DA_flag && (length(data.Date)>=V.NObs*0.5)
            t2=tic;
            [par_update,x0_real_temp,DA.logw_real,neff,beta3_reg_old_real]=da_update(V,DA.logw_real,PAR_real,x0_real_temp,LogL_loc_reg_real,LogL_naz_reg_real,iterDA,data,beta3_reg_old_real,NIter);
            disp(['DA update, iteration ',num2str(iterDA),', neff =', num2str(neff),'; computed in ',num2str(toc(t2)),' s'])
            %x0_real=x0_real_temp;
            clear x0_real_temp
            %if iterDA<NIter
            if iterDA==NIter
               PAR_real_old=PAR_real; 
            end
            PAR_real=par_update;
            %end
        else
            disp('no DA')
            if iterDA<NIter
                iterDA=NIter;
                par_update=PAR_real;
            end
        end
    end
    clear t_real LogL_real ind_out
       
    disp(['all iterations of PF computed in ',num2str(toc(t1)),' s'])
    
    DA.PAR_real=PAR_real;
    
    if ismember(time_model(1),[datenum(2020,3:12,1),datenum(2020,12,15),datenum(2020,12,31)])
        time_x0=time_model(1);
        tt=find(time_red_out==time_model(1));
        time_b=time_red_out(tt:end);
        beta_r=beta_red_out(:,tt:end,:);
        save(['results/state_par',datestr(time_model(1),'yyyy-mm-dd')],'V','PAR_real','time_b','beta_r','x0_real','time_x0')
    end
    %% perform forecast with scenarios
        
    if V.forecast_days>0
        count_for=count_for+1;
        %store state variables at V.times(tm) for next temporal step and forecast and DA
        x0_real=squeeze(x_real(:,end,:)); % initial state variables for next iteration
        for_tstart=DA.times(end);
        if tm==length(V.times)
            time_model_f=for_tstart:max(V.time_final,for_tstart+V.forecast_days);%+V.forecast_days);
        else
            time_model_f=for_tstart:(for_tstart+V.forecast_days);%+V.forecast_days);
        end
        i1=find(V.data.Date(2:end)>(time_model_f(1)),1,'first')+1;
        i2=find(V.data.Date(2:end)<=(time_model_f(end)),1,'last')+1;
        %disp(' ')
        %disp(['forecast at time:', datestr(for_tstart)])
        %' ; final t: ', datestr(V.times(end)), ...
        
        data.Date=V.data.Date(i1:i2);
        
        if length(data.Date)>0
            data.prov_Hnew= V.data.prov_Hnew(:,(i1:i2)-1);
            data.reg_Hnew = V.data.reg_Hnew (:,(i1:i2)-1);
            data.reg_Rnew = V.data.reg_Rnew (:,(i1:i2)-1);
            data.reg_Dnew = V.data.reg_Dnew (:,(i1:i2)-1);
        end
        
        beta_red_old=squeeze(beta_red_real(:,end,:));
        
        ind_out=length(time_model_f);
        x_real_for=zeros(V.NSample,ind_out,V.neqs);
        t_real_for=zeros(V.NSample,ind_out);
        Rt_real_for=zeros(V.NSample,ind_out);
        Et_real_for=zeros(V.NSample,ind_out);

        beta_red_real_for=zeros(V.NSample,ind_out,V.n);
        %LogL_naz_prov_real=zeros(V.NSample,1);
        LogL_naz_reg_real_for =zeros(V.NSample,1);
        %LogL_loc_prov_real=zeros(V.NSample,V.n_reg);
        LogL_loc_reg_real_for =zeros(V.NSample,V.n_reg);
        
        t4=tic;
        %for    cont_sample=1:V.NSample
        parfor cont_sample=1:V.NSample
            Vt=struct;
            Vt.n=n;
            Vt.gammaQgammaH=gammaQgammaH;
            Vt.gammaAgammaQ=gammaAgammaQ;
            Vt.mob=mob;
            Vt.tmob=tmob;
            Vt.nPAR_model=nPAR_model;
            Vt.n_reg=n_reg;
            Vt.prov_IDreg=prov_IDreg;
            Vt.DA_par_flag=DA_par_flag;
            Vt.tbeta=tbeta;
            Vt.beta_string=beta_string;
            Vt.p=Vp;
            Vt.q=Vq;
            Vt.fracH=fracH;
            Vt.prov2reg=prov2reg;
            Vt.data_Date1=data_Date1;
            Vt.time_first_DA=time_first_DA;
            Vt.epsilonP=epsilonP;
            Vt.Wpi=VWpi; % pseudoinverse of W
            Vt.W=VW;
            Vt.LK=LK;
            Vt.mob_red_w=mob_red_w;
            Vt.mob_red_d=mob_red_d;
  
      
            % select current parameters and state variables for the realization
            PAR=PAR_real_old(cont_sample,:);
            Vt.x0=x0_real(cont_sample,:);
            Vt.ProvBetaInc=1;
            
            Vt.beta3_reg_old=beta_red_old(cont_sample,:)';
            Vt.time_model=time_model_f;
            
            %% run model
            [LogPosterior,x,R0,R_t,e_t,beta_red]=sepia_spatial_DA(PAR,Vt,data);
            
            %RR=x(1:i1,8*V.n+1:9*V.n)'; %recorder recovered
            %D=x(1:i1,9*V.n+1:10*V.n)';
            %prov_cumH=x(1:i1,10*V.n+1:11*V.n)'; %cumulative cases
            
            %ind_out=length(Vt.time_model);
            
            x_real_for(cont_sample,:,:)=x;
            t_real_for(cont_sample,:)  =Vt.time_model;
            
            Rt_real_for(cont_sample,:) =R_t; 
            Et_real_for(cont_sample,:) =e_t;
    
            beta_red_real_for(cont_sample,:,:)=beta_red;
            %V.LogL_real.prov(cont_sample,count:(count+ind_out-1))=LogPosterior.prov*ones(1,ind_out);
            %LogL_naz_prov_real(cont_sample)  =LogPosterior.prov;
            LogL_naz_reg_real_for (cont_sample)  =LogPosterior.reg;
            %LogL_loc_prov_real(cont_sample,:)=LogPosterior.loc_prov;
            LogL_loc_reg_real_for (cont_sample,:)=LogPosterior.loc_reg;
        end
        
        disp(['forecast computed in ',num2str(toc(t4)),' s'])
        
        %ind_out=size(x_real,2);
        
        %% extract variables to plot and compute quantiles
        DaTimes=mean(t_real_for,1);
        
        [regHnew_quantile,regR_quantile,regSfrac_quantile,countryHnew_quantile,countryR_quantile,countryE_quantile,regBeta_quantile]=extract_quantiles(PAR_real,V,x_real_for,DaTimes,beta_red_real_for,Rt_real_for,Et_real_for);

        %DA.x_real(:,(count+1):(count+ind_out),:)=x_real;
        DAfor.regHnew_q=regHnew_quantile;
        
        DAfor.countryHnew_q=countryHnew_quantile;
        DAfor.times = DaTimes;
        DAfor.regR_q=regR_quantile;
        DAfor.regBeta_q=regBeta_quantile;
        
        DAfor.countryR=countryR_quantile;
        DAfor.countryE=countryE_quantile;
        
        DAfor.regSfrac_q=regSfrac_quantile;
        DAfor.LogL_naz_reg_real  =LogL_naz_reg_real;
        DAfor.LogL_loc_reg_real  =LogL_loc_reg_real;
        
        
        %% plot forecast scenarios.
        if V.show_fig
            V.nfig_forc=V.nfig_forc+1;
            V.nfig=V.nfig+1;
            plot_DA_results_forc(DA,DAfor,V,time_model)
        else
            if tm==length(V.times)
                plot_DA_results_forc(DA,DAfor,V,time_model)
            end
        end
        %disp('plot forecast finished')
        L_naz_for1=[L_naz_for1,LogL_naz_reg_real];
        L_reg_for1=[L_reg_for1,LogL_loc_reg_real];
    end
    
    %store state variables at V.times(tm) for next temporal step and forecast and DA
    
    if tm>1
        x0_real=squeeze(x_real(:,2,:)); % initial state variables for next iteration
        count=count+1;
    else
        ind=find(DA.times==V.times(1));
        x0_real=squeeze(x_real(:,ind,:));
        count=ind;
%        x0_real=squeeze(x_real(:,end,:));
%       count=length(DA.times);
    end
    %% save/print forecast
    %if V.forecast_days>0
    %    save([folder_name,'/for',datestr(V.times(tm))],'forecast_x_real',...
    %        'forecast_t_real')
    %end
   
end
L_reg1=DA.LogL_loc_reg_real;
L_naz1=DA.LogL_naz_reg_real;

disp(' ')
disp(['all time iterations computed in ',num2str(toc(t0)),' s'])
disp('----------------')
disp(' ')
save([V.folder_name,'/DA_output_',V.filename],'DA','V')
    
%% plot ensemble results with measurement errors
%%
%disp('final plots')

%plot_DA_results(DA,V)
save(['results/beta_red_out',datestr(time_model(end),'yyyy-mm-dd')],'time_red_out','beta_red_out')

return
end