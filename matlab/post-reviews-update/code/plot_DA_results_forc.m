
function plot_DA_results_forc(DA,DAfor,V,time_model)

set(0,'DefaultFigureWindowStyle','docked') %figures are docked
%disp(' in plot_DA_results_forc')
RegFig1=[3 8 1 5 11];
RegFig2=[2 4 6 7 9 10 12:20];

indcol=[1,5,2];



%% extract variables to plot and compute quantiles
%DaTimes=DA.t_real(1,:);
%x_real=DA.x_real;
%beta_red_real=DA.beta_red_real;
%[regHnew_quantile,regR_quantile,regSfrac_quantile,countryHnew_quantile]=extract_quantiles(DA,V,x_real,DaTimes,beta_red_real);

DaTimes_for=DAfor.times;
DaTimes=DA.times;

regHnew_quantile=DA.regHnew_q;
countryHnew_quantile=DA.countryHnew_q;

regHnew_quantile_for=[regHnew_quantile(:,end,:), DAfor.regHnew_q];
countryHnew_quantile_for=[countryHnew_quantile(end,:); DAfor.countryHnew_q];

regR_quantile=DA.regR_q;
regR_quantile_for= DAfor.regR_q;

regB_quantile=DA.regBeta_q;
regB_quantile_for= DAfor.regBeta_q;

%vxlim=[DaTimes(1) max(DaTimes_for(end),V.data.Date(end))];

vxlim=[DaTimes(end)-40,DaTimes_for(end)+5];


%if sum(V.RegSelected)==20
%    RegSuffix='All';
%elseif sum(V.RegSelected)==1
%    RegSuffix=num2str(find(V.RegSelected));
%else
%    error
%end
%filename=['scenario_',V.PostSuffix,'_Reg_',RegSuffix,'_BetaInc_',num2str(V.betaInc*100),'_mob2_',num2str(V.mob_phase2*100)];
%filename=[V.folder_name,'scenario_forc'];

%save(filename,'countryHnew_quantile','regR_quantile','regSfrac_quantile','regHnew_quantile','V','filename','DaTimes')

%%
%figure('Name',['Hospitalized'],'color','w','Visible','Off')
%figure('Name',['Hospitalized'],'color','w')
fig1=figure(4*mod(V.nfig-1,20)+1);
if V.nfig>=20
    clf(fig1)
    fig1=figure(4*mod(V.nfig-1,20)+1);
end

hold OFF
set(fig1,'color','white')% [0 0 1000 900],'color','white')
% results at national level
sb=subplot(3,2,1);
% plot past results
patch([DaTimes(1:end) flip(DaTimes(1:end))]',[countryHnew_quantile(:,1); flip(countryHnew_quantile(:,5))],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
hold on
patch([DaTimes(1:end) flip(DaTimes(1:end))]',[countryHnew_quantile(:,2); flip(countryHnew_quantile(:,4))],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
plot(DaTimes(1:end),countryHnew_quantile(:,3),'color',V.col(1,:),'linewidth',1.5)

% plot forecast
patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',[squeeze(countryHnew_quantile_for(:,1)); flip(squeeze(countryHnew_quantile_for(:,5)))],V.col(3,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
hold on
patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',[squeeze(countryHnew_quantile_for(:,2)); flip(squeeze(countryHnew_quantile_for(:,4)))],V.col(3,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')

plot(DaTimes_for(1:end),squeeze(countryHnew_quantile_for(:,3)),'color',V.col(3,:),'linewidth',1.5)

% plot data
plot(V.data.Date(2:end),sum(V.data.reg_Hnew),'o','color',V.col(2,:),'markersize',4,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')



% plot last DA window
%vylim=[0,max(sum(V.data.reg_Hnew))*1.5];
vylim=[0,max(sum(V.data.reg_Hnew(:,ismember(V.data.Date(2:end),[vxlim(1):vxlim(end)]))))*1.5];
plot([time_model(1), time_model(1)], vylim*0.8,'--k');
plot([time_model(end), time_model(end)], vylim*0.8,'--k');

for ii=3:6
plot([V.tmob(ii), V.tmob(ii)], vylim*0.8,'-g');
end

xticks(sort([datenum(2020,1:12,1),datenum(2020,1:12,15)]))

set(gca,'Xlim',vxlim,'Ylim',vylim)
datetick('x','dd-mm','keeplimits','keepticks')
ylabel('New Hospitalized Cases')
title('Italy')

% pos = get(sb, 'position');
% dim = [pos(1),pos(2)+0.8*pos(4), pos(3)*0.5 0.2*pos(4)];
% str = sprintf('maxLogL =%5.2f\nmeanLogL=%5.2f',max(DAfor.LogL_naz_reg_real),mean(DAfor.LogL_naz_reg_real));
% annotation( 'textbox', dim, 'String', str, 'FitBoxToText','on',  'verticalalignment', 'bottom');


for i=1:length(RegFig1)
    % results at regional level
    r=RegFig1(i);
    sb=subplot(3,2,i+1);
    % past results
    patch([DaTimes(1:end) flip(DaTimes(1:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes(1:end) flip(DaTimes(1:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes(1:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    
    %plot forecast
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regHnew_quantile_for(r,:,1))'; flip(squeeze(regHnew_quantile_for(r,:,5)))'],V.col(3,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regHnew_quantile_for(r,:,2))'; flip(squeeze(regHnew_quantile_for(r,:,4)))'],V.col(3,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes_for(1:end),squeeze(regHnew_quantile_for(r,:,3)),                   'color',V.col(3,:),'linewidth',1.5)
    
    % plot data
    plot(V.data.Date(2:end),V.data.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',4,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    
    vylim=[0,max(max(V.data.reg_Hnew(r,ismember(V.data.Date(2:end),[vxlim(1):vxlim(end)])))*1.5,20)];
    plot([time_model(1), time_model(1)], vylim*0.8,'--k');
    plot([time_model(end), time_model(end)], vylim*0.8,'--k');

    for ii=3:6
        plot([V.tmob(ii), V.tmob(ii)], vylim*0.8,'-g');
    end

    grid on
    set(gca,'Xlim',vxlim,'Ylim',vylim)
    datetick('x','dd-mm','keeplimits','keepticks')

    
    title(V.reg_name(r))
    
%     pos = get(sb, 'position');
%     dim = [pos(1),pos(2)+0.8*pos(4), pos(3)*0.5 0.2*pos(4)];
%     str = sprintf('maxLogL =%5.2f\nmeanLogL=%5.2f',max(DAfor.LogL_loc_reg_real(:,r)),mean(DAfor.LogL_loc_reg_real(:,r)));
%     annotation( 'textbox', dim, 'String', str, 'FitBoxToText','on',  'verticalalignment', 'bottom');

end
if V.save_fig
   % print(['fig/fig1_',filename],'-depsc','-painters')
    saveas(fig1,['fig/Hosp1_',V.filename,'_num',num2str(V.nfig,'%.4d'),'.jpg'],'jpeg')

   %print(['fig/fig1_',V.filename],'-dpng','-r300')
end
%%
%figure('Name',['Hospitalized - others'],'color','w','Visible','Off')
%figure('Name',['Hospitalized - others'],'color','w')

fig2=figure(4*mod(V.nfig-1,20)+2);
if V.nfig>=20
    clf(fig2)
    fig2=figure(4*mod(V.nfig-1,20)+2);
end

hold OFF
set(fig2,'color','w')%'Position',[0 0 1000 1200])
for i=1:20%length(RegFig2)
    % results at regional level
    r=i;%RegFig2(i);
    sb=subplot(5,4,i);
    % past results
    patch([DaTimes(1:end) flip(DaTimes(1:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes(1:end) flip(DaTimes(1:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes(1:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    
    %plot forecast
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regHnew_quantile_for(r,:,1))'; flip(squeeze(regHnew_quantile_for(r,:,5)))'],V.col(3,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regHnew_quantile_for(r,:,2))'; flip(squeeze(regHnew_quantile_for(r,:,4)))'],V.col(3,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes_for(1:end),squeeze(regHnew_quantile_for(r,:,3)),                   'color',V.col(3,:),'linewidth',1.5)
    
    % plot data
    plot(V.data.Date(2:end),V.data.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',4,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    
    vylim=[0,max(max(V.data.reg_Hnew(r,ismember(V.data.Date(2:end),[vxlim(1):vxlim(end)])))*1.5,20)];
    plot([time_model(1), time_model(1)], vylim*0.8,'--k');
    plot([time_model(end), time_model(end)], vylim*0.8,'--k');

    for ii=3:6
        plot([V.tmob(ii), V.tmob(ii)], vylim*0.8,'-g');
    end
    
    grid on
 
    xticks(sort([datenum(2020,1:12,1),datenum(2020,1:12,15)]))
    
    set(gca,'Xlim',vxlim,'Ylim',vylim)
    
    datetick('x','dd-mm','keeplimits','keepticks')
    box off
    title(V.reg_name(r))
    
%     pos = get(sb, 'position');
%     dim = [pos(1),pos(2)+0.6*pos(4), pos(3)*0.5 0.2*pos(4)];
%     str = sprintf('maxLogL =%5.2f\nmeanLogL=%5.2f',max(DAfor.LogL_loc_reg_real(:,r)),mean(DAfor.LogL_loc_reg_real(:,r)));
%     annotation( 'textbox', dim, 'String', str, 'FitBoxToText','on',  'verticalalignment', 'bottom');

end
if V.save_fig
    %print(['fig/fig2_',filename],'-depsc','-painters')
    %print(['fig/fig2_',V.filename],'-dpng','-r300')
    saveas(fig2,['fig/Hosp2_',V.filename,'_num',num2str(V.nfig,'%.4d'),'.jpg'],'jpeg')
end

%% editing
% 
% %plot Sfrac
% 
% 
% figure('Name','S fraction','color','white','Units','Normalized','Position',[0 0 0.5 1])
% for r=1:V.n_reg
%     
%     subplot(5,4,r)
%     patch([DaTimes flip(DaTimes)]',...
%         [squeeze(regSfrac_quantile(r,:,1))'; flip(squeeze(regSfrac_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
%     hold on
%     patch([DaTimes flip(DaTimes)]',...
%         [squeeze(regSfrac_quantile(r,:,2))'; flip(squeeze(regSfrac_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
%     plot(DaTimes,squeeze(regSfrac_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
%     
%     
%     set(gca,'Xlim',vxlim,'Ylim',[0.80 1],'Ygrid','on')
%     if not(mod(r,4)==1)
%         set(gca,'YTicklabel',[])
%     end
%     if r<17
%         set(gca,'XTicklabel',[])
%     else
%         datetick('x','dd-mm','keeplimits','keepticks')
%     end
%     title(V.reg_name(r))
%     box on
% end
% if V.save_fig
%     print(['fig/fig_S_',filename],'-dpng','-r300')
% end
% 
% 
% %% editing

%plot Rt


%figure('Name','R(t)','color','white','Units','Normalized','Position',[0 0 0.5 1])
fig3=figure(4*mod(V.nfig-1,20)+3);
if V.nfig>=20
    clf(fig3)
    fig3=figure(4*mod(V.nfig-1,20)+3);
end
set(fig3,'color','white')

yl=max(max(regR_quantile(:,:,5)));
yl=max(yl,1.2);
vylim=[0,3];
for r=1:V.n_reg
    
    subplot(5,4,r)
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,1))'; flip(squeeze(regR_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,2))'; flip(squeeze(regR_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes,squeeze(regR_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    plot([DaTimes(1) DaTimes(end)+30],[1 1],'--r')
    
    grid on
      %plot forecast

    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regR_quantile_for(r,:,1))'; flip(squeeze(regR_quantile_for(r,:,5)))'],V.col(3,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regR_quantile_for(r,:,2))'; flip(squeeze(regR_quantile_for(r,:,4)))'],V.col(3,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes_for(1:end),squeeze(regR_quantile_for(r,:,3)),                   'color',V.col(3,:),'linewidth',1.5)
    set(gca,'Xlim',vxlim,'Ylim',vylim)
   plot([time_model(1), time_model(1)], vylim*0.8,'--k');
    plot([time_model(end), time_model(end)], vylim*0.8,'--k');

   
    if not(mod(r,4)==1)
        set(gca,'YTicklabel',[])
    elseif r==9
        ylabel('estimate of regional R_t')
    end
    if r<17
        set(gca,'XTicklabel',[])
    else
        datetick('x','dd-mm','keeplimits','keepticks')
    end
    title(V.reg_name(r))
    %box on
    
    yyaxis right
    plot(V.mob_red_d,V.mob_red_w_reg(:,r)/100)    
   ylim([-0.8,0.1])
   yticks([-0.8 -0.5 -0.2 1.0'])
    if not(mod(r,4)==0)
        set(gca,'YTicklabel',[])
    else
        if r==12
        ylabel('google mobility reduction')
        end
    end
    
 
    
end

%%
%if V.save_fig
 saveas(fig3,['fig/R_',V.filename,'_num',num2str(V.nfig,'%.4d'),'.jpg'],'jpeg')

    %print(['fig/R_',V.filename],'-dpng','-r300')
%end
fig4=figure(4*mod(V.nfig-1,20)+4);
if V.nfig>=20
    clf(fig4)
    fig4=figure(4*mod(V.nfig-1,20)+4);
end
set(fig4,'color','white')

yl=max(max(regB_quantile(:,:,5)));
%yl=max(yl,1.2);
vylim=[0,1.2];
for r=1:V.n_reg
    
    subplot(5,4,r)
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regB_quantile(r,:,1))'; flip(squeeze(regB_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regB_quantile(r,:,2))'; flip(squeeze(regB_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes,squeeze(regB_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
   % plot([DaTimes(1) DaTimes(end)+30],[1 1],'--r')
    
    grid on
      %plot forecast

    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regB_quantile_for(r,:,1))'; flip(squeeze(regB_quantile_for(r,:,5)))'],V.col(3,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes_for(1:end) flip(DaTimes_for(1:end))]',...
        [squeeze(regB_quantile_for(r,:,2))'; flip(squeeze(regB_quantile_for(r,:,4)))'],V.col(3,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes_for(1:end),squeeze(regB_quantile_for(r,:,3)),                   'color',V.col(3,:),'linewidth',1.5)
    set(gca,'Xlim',vxlim,'Ylim',vylim)
        yticks([0 0.4 0.8 1.2])

   plot([time_model(1), time_model(1)], vylim*0.8,'--k');
    plot([time_model(end), time_model(end)], vylim*0.8,'--k');

   
    if not(mod(r,4)==1)
        set(gca,'YTicklabel',[])
    elseif r==9
        ylabel('coefficient of \beta_0')
    end
    if r<17
        set(gca,'XTicklabel',[])
    else
        datetick('x','dd-mm','keeplimits','keepticks')
    end
    title(V.reg_name(r))
    %box on
    
    yyaxis right
    plot(V.mob_red_d,V.mob_red_w_reg(:,r)/100)    
   ylim([-0.8,0.1])
    yticks([-0.8 -0.5 -0.2 1.0'])
   
    if not(mod(r,4)==0)
        set(gca,'YTicklabel',[])
    else
        if r==12
        ylabel('google mobility reduction')
        end
    end    
end
 saveas(fig4,['fig/Beta_',V.filename,'_num',num2str(V.nfig,'%.4d'),'.jpg'],'jpeg')


return 