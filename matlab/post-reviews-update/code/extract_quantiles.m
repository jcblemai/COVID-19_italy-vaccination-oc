%% extract variables to plot and compute quantiles
function [regHnew_quantile,regR_quantile,regSfrac_quantile,countryHnew_quantile,countryR_quantile,countryE_quantile,regBeta_quantile]=extract_quantiles(PAR_real,V,x_real,DaTimes,beta_red_real,Rt_real,Et_real)
cumH_real=x_real(:,:,10*V.n+1:11*V.n);
prov_S_real=x_real(:,:,1:V.n);
dt=diff(DaTimes);
ind1=find(dt==0);

nt=size(x_real,2);
reg_R_real=zeros(V.n_reg,nt,V.NSample);
reg_Sfrac_real=zeros(V.n_reg,nt,V.NSample);
reg_beta_red_real=zeros(V.n_reg,nt,V.NSample);

reg_Hnew_real=zeros(V.n_reg,nt-1,V.NNoise4sample,V.NSample);
N=V.N;
n=V.n;
prov2reg=V.prov2reg;
prov_ave_reg=V.prov_ave_reg;
NNoise4sample=V.NNoise4sample;
R0_real=V.R0_real;
parfor cont_sample=1:V.NSample
%for cont_sample=1:V.NSample
    omega=PAR_real(cont_sample,12);
    
    cumH=squeeze(cumH_real(cont_sample,:,:)); %cumulative cases
    dsumH=diff(cumH);
    
    prov_Hnew=dsumH./repmat(dt',1,n);
    prov_Hnew(ind1,:)=prov_Hnew(ind1-1,:);
    
    prov_S=squeeze(prov_S_real(cont_sample,:,:))';
    
    reg_Sfrac=prov2reg*prov_S./(prov2reg*N);
    
    reg_R_real(:,:,cont_sample)=R0_real(cont_sample)*(prov_ave_reg*(squeeze(beta_red_real(cont_sample,:,:))')).*reg_Sfrac;
    reg_beta_red_real(:,:,cont_sample)=prov_ave_reg*(squeeze(beta_red_real(cont_sample,:,:))');
    
    reg_Sfrac_real(:,:,cont_sample)=reg_Sfrac;
    
    r=(prov_Hnew')./(omega-1);
    p=1./omega;
    reg_Hnew_real1=[];
    for cont_noise=1:NNoise4sample
        %cont_real=cont_real+1;
        %add error at the prov_Hnew_sim
        reg_Hnew_real1(:,:,cont_noise)=prov2reg*nbinrnd(r,p);
    end
    reg_Hnew_real(:,:,:,cont_sample)=reg_Hnew_real1;
    %reg_Hnew_real1(:,:,cont_real)=V.prov2reg*prov_Hnew';
end

reg_Hnew_real=reshape(reg_Hnew_real,V.n_reg,nt-1,V.NSample*V.NNoise4sample);

regHnew_quantile=quantile(reg_Hnew_real,V.quant,3);
regR_quantile=quantile(reg_R_real,V.quant,3);
regBeta_quantile=quantile(reg_beta_red_real,V.quant,3);
regSfrac_quantile=quantile(reg_Sfrac_real,V.quant,3);
countryHnew_quantile=quantile(squeeze(sum(reg_Hnew_real)),V.quant,2);
countryR_quantile=quantile(Rt_real(:,2:end),V.quant,1)';
countryE_quantile=quantile(Et_real(:,2:end),V.quant,1)';

%%
return
end