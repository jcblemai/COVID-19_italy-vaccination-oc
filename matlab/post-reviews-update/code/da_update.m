function [par,statex,logw,neff,beta3_reg_old_real]=da_update(V,LogWIn,PAR_real,statex,...
    LogL_loc_reg_real,LogL_naz_reg_real,iterDA,data,beta3_reg_old_real,NIter)
                         
par=PAR_real;
loglike=LogL_naz_reg_real;
logw=LogWIn; %log weights
NSample=V.NSample;
resampling_loc=1;
[w]=update_weights(logw,loglike);

% check if resampling is needed
neff=1/sum(w.^2);
%disp(['% neff=', num2str(neff/NSample)])


if neff<V.rtau*NSample
    %disp(['systematic resampling'])
    
    %
    % decrease std_beta3 at each DA iteration.
    % in the last iteration, perturbation use std_beta3 for the next iteration.
    if iterDA<NIter
        std_beta3=(V.std_beta3f-V.std_beta3i)/(NIter-1) *(iterDA)+V.std_beta3i;
        disp(std_beta3)
    else
        std_beta3=V.std_beta3i;
    end
    %disp(std_beta3)
    %     if resampling_loc==0
    %         [indn]=systematic_resampling(w,NSample);
    %         logw=-ones(NSample,1)*log(NSample);
    %
    %         % update state variables
    %         oldx=statex;
    %         statex=oldx(indn,:);
    %         % the following allows that resampled particles stay associated with their
    %         % parameters (not needed if states are not updated).
    %
    %         [resampledI,IA,JA]=unique(indn,'first');
    %         %     code for testing
    %         %     count=1;
    %         %     v=1:NSample;
    %         %     indo=v;
    %         %     v(IA)=indo(resampledI);
    %         %
    %         % keep parameters coupled with the updated variables
    %         parold=par;
    %         par(IA,:)=parold(resampledI,:);
    %
    %         count=1;
    %         for i=1:length(resampledI)
    %             if ~ismember(resampledI(i),IA)
    %                 while ismember(IA(count),(resampledI))
    %                     count=count+1;
    %                 end
    %                 par(resampledI(i),:)=parold(IA(count),:);
    %                 %v(resampledI(i))=indo(IA(count));
    %                 count=count+1;
    %             end
    %         end
    %
    %         if V.DA_par_flag
    %             % update beta3
    %
    %             for nr=1:V.n_reg
    %                 % perturb parameters of particles that were not duplicated
    %                 for cont_sample=1:NSample
    %                     if ~ismember(cont_sample,resampledI)
    %                         mean_beta3=par(indn(cont_sample),V.nPAR_model+nr);
    %                         %gt=makedist('Normal','mu',mean_beta3,'sigma',std_beta3);
    %                         %gt=truncate(gt,0,1.1);
    %
    %
    %                         gt=makedist('Normal','mu',1,'sigma',std_beta3);
    %                         gt=truncate(gt,0.8,1.2);
    %                         beta3=mean_beta3*random(gt,1,1);
    %                         par(cont_sample,V.nPAR_model+nr)=beta3;
    %                     end
    %                 end
    %             end
    %         end
    %     else
    
    par1=zeros(NSample,V.n_reg);
    par_temp=par(:,(V.nPAR_model+1):(V.nPAR_model+V.n_reg));
    temp_beta3_old=zeros(V.NSample,V.n_reg);
    for reg=1:V.n_reg
        indp=find(V.prov_IDreg==reg);
        beta3_old_old(:,reg)=beta3_reg_old_real(:,indp(1));
    end
    
    parfor reg=1:V.n_reg
    %for reg=1:V.n_reg
        %perform local update of the parameters
        %if mean(data.reg_Hnew(reg,:))>1
        loglike=LogL_loc_reg_real(:,reg);
        logw=LogWIn; %log weights
        
        [w]=update_weights(logw,loglike);
        [indnpar]=systematic_resampling(w,NSample);
        [resampledI,IA,JA]=unique(indnpar,'first');
        
        indp=find(V.prov_IDreg==reg);
        % perturb parameters of duplicated particles
        beta3=zeros(NSample,1);
        
        for cont_sample=1:NSample
            mean_beta3=min(max(par_temp(indnpar(cont_sample),reg),0.2),1.5);
            %mean_beta3=par_temp(indnpar(cont_sample),reg);
            temp_beta3_old(cont_sample,reg)=beta3_old_old(indnpar(cont_sample),reg);
            if ~ismember(cont_sample,resampledI)
                gt=makedist('Normal','mu',mean_beta3,'sigma',std_beta3);
                gt=truncate(gt,V.min_truncate*mean_beta3,V.max_truncate*mean_beta3);
                
                %gt=makedist('Normal','mu',1,'sigma',std_beta3);
                %gt=truncate(gt,V.min_truncate,V.max_truncate);
                
                beta3(cont_sample)=random(gt,1,1);
                %beta3(cont_sample)=mean_beta3*random(gt,1,1);
            else
                beta3(cont_sample)=mean_beta3;
            end
        end
        par1(:,reg)=beta3;
    end
    
    for reg=1:V.n_reg
        indp=find(V.prov_IDreg==reg);
        beta3_reg_old_real(:,indp)=repmat(temp_beta3_old(:,reg),1,length(indp));
    end
    
    par(:,(V.nPAR_model+1):(V.nPAR_model+V.n_reg))=par1;
    %end
    %     end
    
    logw=-ones(NSample,1)*log(NSample);
else
    disp(['NO RESAMPLING in PF'])
    logw=log(w);
end

return
end
