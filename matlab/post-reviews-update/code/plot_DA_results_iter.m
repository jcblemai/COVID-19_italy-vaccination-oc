function plot_DA_results_iter(DA,V,iterDA,time_model)

set(0,'DefaultFigureWindowStyle','docked') %figures are docked
RegFig1=[3 8 1 5 11];
RegFig2=[2 4 6 7 9 10 12:20];

DaTimes=DA.times;
%% extract variables to plot and compute quantiles

Hnew=[];
cont_real=0;
vxlim=[DaTimes(1) V.data.Date(end)];
vxlim=[DaTimes(end)-50,DaTimes(end)+10];
%for cont_sample=1:V.NSample
    %omega=DA.PAR_real(cont_sample,12);
    
    %cumH=squeeze(DA.x_real(cont_sample,:,10*V.n+1:11*V.n)); %cumulative cases
    %dsumH=diff(cumH);
    %dt=diff(DaTimes);
    %ind1=find(dt==0);
    %prov_Hnew=dsumH./repmat(dt',1,V.n);
    %prov_Hnew(ind1,:)=prov_Hnew(ind1-1,:);
 
 %   prov_S=squeeze(DA.x_real(cont_sample,:,1:V.n))';
    
 %   reg_Sfrac=V.prov2reg*prov_S./(V.prov2reg*V.N);
  %  reg_R_real(:,:,cont_sample)=DA.R0_real(cont_sample)*(V.prov_ave_reg*(squeeze(DA.beta_red_real(cont_sample,:,:))')).*reg_Sfrac;
  %  reg_Sfrac_real(:,:,cont_sample)=reg_Sfrac;
    
    %r=(prov_Hnew')./(omega-1);
    %p=1./omega;
    %for cont_noise=1:DA.NNoise4sample
    %    cont_real=cont_real+1;
    %    %add error at the prov_Hnew_sim
    %    reg_Hnew_real(:,:,cont_real)=V.prov2reg*nbinrnd(r,p);    
    %end
    %reg_Hnew_real1(:,:,cont_real)=V.prov2reg*prov_Hnew';
%end

%regHnew_quantile=quantile(reg_Hnew_real,V.quant,3);
regR_quantile=DA.regR_q;
%regSfrac_quantile=quantile(reg_Sfrac_real,V.quant,3);
%countryHnew_quantile=quantile(squeeze(sum(reg_Hnew_real)),V.quant,2);
%%

%if sum(V.RegSelected)==20
%    RegSuffix='All';
%elseif sum(V.RegSelected)==1
%    RegSuffix=num2str(find(V.RegSelected));
%else
%    error
%end
%filename=['scenario_',V.PostSuffix,'_Reg_',RegSuffix,'_BetaInc_',num2str(V.betaInc*100),'_mob2_',num2str(V.mob_phase2*100)];
filename=V.filename;

%save(['posterior/',filename],'countryHnew_quantile','regR_quantile','regSfrac_quantile','regHnew_quantile','V','filename','DaTimes')

%%
%figure('Name',['Hospitalized'])
%set(gcf,'color','w')

% subplot(3,2,1)
% patch([DaTimes(2:end) flip(DaTimes(2:end))]',[countryHnew_quantile(:,1); flip(countryHnew_quantile(:,5))],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
% hold on
% patch([DaTimes(2:end) flip(DaTimes(2:end))]',[countryHnew_quantile(:,2); flip(countryHnew_quantile(:,4))],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
% plot(V.Date(2:end),sum(V.reg_Hnew),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
% plot(DaTimes(2:end),countryHnew_quantile(:,3),'color',V.col(1,:),'linewidth',1.5)
% box off
% set(gca,'Xlim',[DaTimes(2) DaTimes(end)])
% datetick('x','dd-mm','keeplimits','keepticks')
% ylabel('New Hospitalized Cases')
% title('Italy')
% 
% for i=1:length(RegFig1)
%     r=RegFig1(i);
%     subplot(3,2,i+1)
%     patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
%         [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
%     hold on
%     patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
%         [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
%     plot(V.Date(2:end),V.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
%     plot(DaTimes(2:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
%     box off
%     set(gca,'Xlim',[DaTimes(2) DaTimes(end)])
%     datetick('x','dd-mm','keeplimits','keepticks')
%     box off
%     title(V.reg_name(r))
%     
% end
% if V.save_fig
%    % print(['fig/fig1_',filename],'-depsc','-painters')
%     print(['fig/fig1_',filename],'-dpng','-r300')
% end
%%
%% editing

%plot Sfrac

%plot Rt
%fig=figure('Name',['R(t) estimate, iteration',num2str(iterDA)],'color','white','Units','Normalized','Visible','Off');
%fig=figure('Name',['R(t) estimate, iteration',num2str(iterDA)],'color','white','Units','Normalized');

fig=figure(mod(V.nfig,20)+1);
if V.nfig>=20
    clf(fig)
    fig=figure(mod(V.nfig,20)+1);
end

%clear('fig')
hold OFF
for r=1:V.n_reg    
    sb=subplot(5,4,r);
    if iterDA<V.NIter
        patch([time_model flip(time_model)], [0*time_model 0*time_model+4],[1,0,0],'edgecolor','n','facealpha',0.05,'Handlevisibility','off')
    else
        patch([time_model flip(time_model)], [0*time_model 0*time_model+4],[0,1,0],'edgecolor','n','facealpha',0.05,'Handlevisibility','off')
    end
    hold on
    
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,1))'; flip(squeeze(regR_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,2))'; flip(squeeze(regR_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes,squeeze(regR_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    % plot last DA window
vylim=[0,3.5];

for ii=3:6
plot([V.tmob(ii), V.tmob(ii)], vylim,'-g');
end
    xticks(sort([datenum(2020,1:12,1),datenum(2020,1:12,15)]))
    yticks([1,2,3])
    set(gca,'Xlim',vxlim,'Ylim',vylim,'Ygrid','on')
    if not(mod(r,4)==1)
        set(gca,'YTicklabel',[])
    end
    if r<17
        set(gca,'XTicklabel',[])
    else
        
        datetick('x','dd-mm','keeplimits','keepticks')
    end
    %title([V.reg_name(r),',logL=', num2str(max(LogL_real.loc_reg(:,r)),'%5.2f')])
    title(V.reg_name{r});
    box on
    
    
%     pos = get(sb, 'position');
%     dim = pos.*[1 1 0.5 0.5];
%     str = sprintf('maxLogL =%5.2f\nmeanLogL=%5.2f',max(DA.LogL_loc_reg_real(:,end,r)),mean(DA.LogL_loc_reg_real(:,end,r)));
%     annotation( 'textbox', dim, 'String', str, 'FitBoxToText','on',  'verticalalignment', 'bottom');
end

%%
%%
if V.save_fig
    if iterDA==V.NIter
    %print(['fig/Rt_',V.filename],'-dpng','-r300')
    saveas(fig,['fig/Rt_',V.filename,'_num',num2str(V.nfig_iter,'%.4d'),'.jpg'],'jpeg')
    end    
end


return 