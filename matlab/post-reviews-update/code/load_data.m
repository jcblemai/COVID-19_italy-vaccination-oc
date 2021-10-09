function V=load_data(V)

V = rmfield(V,{'Date','prov_C','reg_C','reg_H','reg_R','reg_D','reg_Q',...
    'reg_Hnew','prov_Cnew','reg_Rnew','reg_Dnew','reg_Cnew','reg_Cnew_unknown'});


%load case data
load ../dati_covid_italy/data-model-matlab/data_ProvReg;


V.data.Date=date_vec; clear date_vec
V.data.prov_C=provinceCases(not(isnan(provinceCases(:,1))),:);  %, excluide non exixting ID

% V.prov_Inc=diff(V.prov_C,1,2);  
V.n_reg=20;
V.data.reg_C=regionCases.TOT;  %total cumulative cases by region
V.data.reg_H=regionCases.INF+regionCases.ICU;
V.data.reg_R=regionCases.REC;
V.data.reg_D=regionCases.DEAD;
V.data.reg_Q=regionCases.HOME;
V.data.reg_Hnew=regionCases.NEWHOSP(:,2:end);


V.data.prov_Cnew=diff(V.data.prov_C,1,2);
V.data.reg_Rnew=diff(V.data.reg_R,1,2);
V.data.reg_Dnew=diff(V.data.reg_D,1,2);
V.data.reg_Cnew=diff(V.data.reg_C,1,2);

%check H cases do not exceed 
%V.data.reg_Hnew=min(V.data.reg_Hnew,V.data.reg_Cnew);

V.data.reg_Cnew_unknown=V.prov2reg*V.data.prov_Cnew;  %I re_calculate this this way so the downscaling using prow_Cnew is ok (prov_Cnew do not have the repatition of the unkwown cases)

% Guess new hospitalized at province level
V.data.prov_Hnew=zeros(size(V.data.prov_Cnew));

for t=1:size(V.data.prov_Cnew,2)
    for i=1:V.n
        if V.data.reg_Hnew(V.prov_IDreg(i),t)>0 && V.data.reg_Cnew_unknown(V.prov_IDreg(i),t)>0
            V.data.prov_Hnew(i,t)=round(V.data.reg_Hnew(V.prov_IDreg(i),t)*V.data.prov_Cnew(i,t)/V.data.reg_Cnew_unknown(V.prov_IDreg(i),t));
        end
    end
end


if sum(sum(isnan(V.data.prov_Hnew)))>0
    error('NaN in prov_Hnew');
end
if sum(sum((V.data.prov_Hnew<0)))>0
     error('negative number in prov_Hnew');
end

V.data.prov_Hnew_nosmooth=V.data.prov_Hnew;
V.data.prov_Hnew=movmean(V.data.prov_Hnew,7,2);

clear t regionCases regionCases_raw provinceCases ...
    provinceCases_raw par i