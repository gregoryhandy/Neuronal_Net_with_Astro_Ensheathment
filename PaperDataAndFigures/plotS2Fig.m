%%
% Plots spatial correlation curves across different values of probability 
% and strength of ensheathment 
%
% Note: only excitatory neurons are ensheathed 
%%


%% load the data and get the statistics

load('./CompressedData/Spatial/Spatial_SuppData_condensed.mat')

%%
figure(9); clf;
set(gca,'fontsize',16)
% temp_title = sprintf('\\sigma_{rec} = %.2f, Prob_{Ensheath} = %.2f',param_matrix(1,1),param_matrix(2,1));

subplot(1,3,1); hold on;
set(gca,'fontsize',16)
temp_title = sprintf('p_{en} = %.2f',param_matrix(2,13));
title(temp_title)
count =1;
vis_vec = [1:-0.9/3:.1];
for jj = 13:16
    hold on
    plot(distbins(jj,:), SCcorrSimMean(jj,:),'linewidth',1.5,'color',[1 0.5 0 vis_vec(count)])
    count = count + 1;
end
x_temp = [0:0.01:1];
plot(x_temp,0*x_temp,'k--')
xlim([0 0.7])
ylim([-0.04 .105])
xlabel('Distance (a.u.)')
ylabel('Spike count correlations')


subplot(1,3,2); hold on;
set(gca,'fontsize',16)
temp_title = sprintf('p_{en} = %.2f',param_matrix(2,17));
title(temp_title)
count =1;
for jj = 17:20
    plot(distbins(jj,:), SCcorrSimMean(jj,:),'linewidth',1.5,'color',[1 0.5 0 vis_vec(count)])
    count = count + 1;
end
x_temp = [0:0.01:1];
plot(x_temp,0*x_temp,'k--')
xlim([0 0.7])
ylim([-0.04 .105])
xlabel('Distance (a.u.)')

subplot(1,3,3); hold on
set(gca,'fontsize',16)
temp_title = sprintf('p_{en} = %.2f',param_matrix(2,21));
title(temp_title)
legend_names_2 = [];
count = 1;
for jj = 21:24
    plot(distbins(jj,:), SCcorrSimMean(jj,:),'linewidth',1.5,'color',[1 0.5 0 vis_vec(count)])
    count = count + 1;
    legend_names_2 = [legend_names_2, sprintf("s_{en}=%.2f",1-param_matrix(3,jj))];
end

x_temp = [0:0.01:1];
plot(x_temp,0*x_temp,'k--')
xlim([0 0.7])
ylim([-0.04 .105])
legend(legend_names_2)
xlabel('Distance (a.u.)')