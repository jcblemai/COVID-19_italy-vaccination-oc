function plot_DA_results(DA,V)

set(0,'DefaultFigureWindowStyle','docked') %figures are docked
RegFig1=[3 8 1 5 11];
RegFig2=[2 4 6 7 9 10 12:20];

DaTimes=DA.t_real(1,:);
%% extract variables to plot and compute quantiles

Hnew=[];
cont_real=0;
for cont_sample=1:V.NSample
    omega=DA.PAR_real(cont_sample,12);
    
    cumH=squeeze(DA.x_real(cont_sample,:,10*V.n+1:11*V.n)); %cumulative cases
    dsumH=diff(cumH);
    dt=diff(DaTimes);
    ind1=find(dt==0);
    prov_Hnew=dsumH./repmat(dt',1,V.n);
    prov_Hnew(ind1,:)=prov_Hnew(ind1-1,:);
 
    prov_S=squeeze(DA.x_real(cont_sample,:,1:V.n))';
    
    reg_Sfrac=V.prov2reg*prov_S./(V.prov2reg*V.N);
    reg_R_real(:,:,cont_sample)=DA.R0_real(cont_sample)*(V.prov_ave_reg*(squeeze(DA.beta_red_real(cont_sample,:,:))')).*reg_Sfrac;
    reg_Sfrac_real(:,:,cont_sample)=reg_Sfrac;
    
    
    r=(prov_Hnew')./(omega-1);
    p=1./omega;
    for cont_noise=1:V.NNoise4sample
        cont_real=cont_real+1;
        %add error at the prov_Hnew_sim
        reg_Hnew_real(:,:,cont_real)=V.prov2reg*nbinrnd(r,p);
        
    end
    %reg_Hnew_real1(:,:,cont_real)=V.prov2reg*prov_Hnew';
end


regHnew_quantile=quantile(reg_Hnew_real,V.quant,3);
regR_quantile=quantile(reg_R_real,V.quant,3);
regSfrac_quantile=quantile(reg_Sfrac_real,V.quant,3);
countryHnew_quantile=quantile(squeeze(sum(reg_Hnew_real)),V.quant,2);
%%

if sum(V.RegSelected)==20
    RegSuffix='All';
elseif sum(V.RegSelected)==1
    RegSuffix=num2str(find(V.RegSelected));
else
    error
end
filename=['scenario_',V.PostSuffix,'_Reg_',RegSuffix,'_BetaInc_',num2str(V.betaInc*100),'_mob2_',num2str(V.mob_phase2*100)];

save(['posterior/',filename],'countryHnew_quantile','regR_quantile','regSfrac_quantile','regHnew_quantile','V','filename','DaTimes')

%%
fig=figure('Name',['Hospitalized']);
set(fig,'color','w','Position',[0 0 1000 900])

subplot(3,2,1)
patch([DaTimes(2:end) flip(DaTimes(2:end))]',[countryHnew_quantile(:,1); flip(countryHnew_quantile(:,5))],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
hold on
patch([DaTimes(2:end) flip(DaTimes(2:end))]',[countryHnew_quantile(:,2); flip(countryHnew_quantile(:,4))],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
plot(V.Date(2:end),sum(V.reg_Hnew),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
plot(DaTimes(2:end),countryHnew_quantile(:,3),'color',V.col(1,:),'linewidth',1.5)
box off
set(gca,'Xlim',[DaTimes(2) DaTimes(end)])
datetick('x','dd-mm','keeplimits','keepticks')
ylabel('New Hospitalized Cases')
title('Italy')

for i=1:length(RegFig1)
    r=RegFig1(i);
    subplot(3,2,i+1)
    patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(V.Date(2:end),V.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    plot(DaTimes(2:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    set(gca,'Xlim',[DaTimes(2) DaTimes(end)])
    datetick('x','dd-mm','keeplimits','keepticks')
    box off
    title(V.reg_name(r))
    
end
if V.save_fig
   % print(['fig/fig1_',filename],'-depsc','-painters')
    print(['fig/fig1_',filename],'-dpng','-r300')
end
%%
fig=figure('Name',['Hospitalized - others']);
set(fig,'color','w','Position',[0 0 1000 1200])
for i=1:length(RegFig2)
    r=RegFig2(i);
    subplot(5,3,i)
    patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
        [squeeze(regHnew_quantile(r,:,1))'; flip(squeeze(regHnew_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes(2:end) flip(DaTimes(2:end))]',...
        [squeeze(regHnew_quantile(r,:,2))'; flip(squeeze(regHnew_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(V.Date(2:end),V.reg_Hnew(r,:),'o','color',V.col(2,:),'markersize',5,'markerfacecolor','w','linewidth',1.5,'displayname','New Hospitalized')
    plot(DaTimes(2:end),squeeze(regHnew_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    box off
    set(gca,'Xlim',[DaTimes(2) DaTimes(end)])
    datetick('x','dd-mm','keeplimits','keepticks')
    box off
    title(V.reg_name(r))
    
end
if V.save_fig
    %print(['fig/fig2_',filename],'-depsc','-painters')
    print(['fig/fig2_',filename],'-dpng','-r300')
end

%% editing

%plot Sfrac


figure('Name','S fraction','color','white','Units','Normalized','Position',[0 0 0.5 1])
for r=1:V.n_reg
    
    subplot(5,4,r)
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regSfrac_quantile(r,:,1))'; flip(squeeze(regSfrac_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regSfrac_quantile(r,:,2))'; flip(squeeze(regSfrac_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes,squeeze(regSfrac_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    set(gca,'Xlim',[DaTimes(2) DaTimes(end)],'Ylim',[0.80 1],'Ygrid','on')
    if not(mod(r,4)==1)
        set(gca,'YTicklabel',[])
    end
    if r<17
        set(gca,'XTicklabel',[])
    else
        datetick('x','dd-mm','keeplimits','keepticks')
    end
    title(V.reg_name(r))
    box on
end
%if V.save_fig
    print(['fig/fig_S_',filename],'-dpng','-r300')
%end


%% editing

%plot Rt


figure('Name','R(t)','color','white','Units','Normalized','Position',[0 0 0.5 1])
for r=1:V.n_reg
    
    subplot(5,4,r)
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,1))'; flip(squeeze(regR_quantile(r,:,5)))'],V.col(1,:),'edgecolor','n','facealpha',0.15,'Handlevisibility','off')
    hold on
    patch([DaTimes flip(DaTimes)]',...
        [squeeze(regR_quantile(r,:,2))'; flip(squeeze(regR_quantile(r,:,4)))'],V.col(1,:),'edgecolor','n','facealpha',0.25,'Handlevisibility','off')
    plot(DaTimes,squeeze(regR_quantile(r,:,3)),'color',V.col(1,:),'linewidth',1.5)
    set(gca,'Xlim',[DaTimes(2) DaTimes(end)],'Ylim',[0 2.5],'Ygrid','on')
    if not(mod(r,4)==1)
        set(gca,'YTicklabel',[])
    end
    if r<17
        set(gca,'XTicklabel',[])
    else
        datetick('x','dd-mm','keeplimits','keepticks')
    end
    title(V.reg_name(r))
    box on
end

%%
if V.save_fig
    print(['fig/R_',filename],'-dpng','-r300')
end


return 