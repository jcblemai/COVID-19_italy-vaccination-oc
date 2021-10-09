function plot_Posterior_sample(Posterior_sample,V)
%
figure('color','w','Position',[0 0 1000 600])
for i=1:V.nPAR_model
    subplot(2,ceil(((V.nPAR_model)/2)),i)
    plot(Posterior_sample(:,i))
    xlabel(V.PAR_names{i})
    box off
end
% 
figure('color','w','Position',[0 0 1000 600])
for i=1:V.nPAR_model
    subplot(3,ceil((V.nPAR_model/3)),i)
    h=histogram(Posterior_sample(:,i),10,'Normalization','pdf','edgecolor','none');
    xlabel(V.PAR_names{i})
    box off
    ylabel('Probability')
end


%%
% figure('color','w','Position',[0 0 1300 800])
% plotmatrix(Posterior_sample(:,1:V.nPAR_model))

end