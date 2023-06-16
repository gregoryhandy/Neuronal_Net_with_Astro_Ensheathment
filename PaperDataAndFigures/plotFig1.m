%%
% Plots the results of DiRT simulations for different levels of 
% astrocyte protrusions 
%%

%% Figure 1B: Panel for the number of activated receptors over time curves

% Load the data
load('./CompressedData/DiRTSimulations/AstroDiRTData.mat')

figure(1); clf; 
color_scheme = [0.85 0.325 0.098; 0 0.447 0.741; 0.466 0.674 0.188; 0 0 0];
ex_indices = [9,13];

subplot(1,3,1); hold on; h=[];
for ii = 1:length(ex_indices)
    
    h(ii) = plot(time_DiRT_all(1:end,ex_indices(ii)),50-r_ave_DiRT_all(1:end,ex_indices(ii)),'linewidth',3,'color',color_scheme(ii,:));
    plot(time_DiRT_all(half_max_start(ex_indices(ii)),ex_indices(ii)), max_val(ex_indices(ii))/2,'.','markersize',20,'color',color_scheme(ii,:))
    plot(time_DiRT_all(half_max_end(ex_indices(ii)),ex_indices(ii)), max_val(ex_indices(ii))/2,'.','markersize',20,'color',color_scheme(ii,:))
    
    x_temp = [time_DiRT_all(half_max_start(ex_indices(ii)))...
        :0.01:time_DiRT_all(half_max_end(ex_indices(ii)))];

    plot(x_temp,max_val(ex_indices(ii))/2*ones(length(x_temp),1),'--','linewidth',2,'color',color_scheme(ii,:))
end
set(gca,'fontsize',16)
xlim([0 1])
legend(h,strcat('p= ',legend_names(ex_indices)))
xlabel('Time (a.u.)')
ylabel('# of Activated Receptors')
yticks([0:10:50])


%% Fit the piecewise linear model for the full range of parameter values

% Fit the model in the two regimes
% Saturated regime
xAUCPart1 = [(1-domain_x_length_vec(1:7)) -0.15];
yAUCPart1 = mean(AUC(1:6))*ones(8,1);

xHMWPart1 = [(1-domain_x_length_vec(1:7)) -0.15];
yHMWPart1 = mean(half_max_width(1:6))*ones(8,1);

% Fit the linear regime for the AUC
[linear_model_AUC,goodness_AUC] = fit((1-domain_x_length_vec(7:end))',AUC(7:end),'poly1');
xAUCPart2 = [-0.15:0.001:1];
fitStart = [xAUCPart1(end), yAUCPart1(end)];
yAUCPart2 = linear_model_AUC.p1*(xAUCPart2-fitStart(1))+fitStart(2);

% Fit the linear regime for the half-max width
[linear_model_HMW,~] = fit((1-domain_x_length_vec(7:end-4))',...
    half_max_width(7:end-4),'poly1');
xHMWPart2 = [-0.15:0.001:1];
fitStart = [xHMWPart1(end), yHMWPart1(end)];
yHMWPart2 = linear_model_HMW.p1*(xHMWPart2-fitStart(1))+fitStart(2);

%% Panels with data and piecewise linear fits
subplot(1,3,2); hold on;
plot((1-domain_x_length_vec),AUC,'ko','linewidth',1.5,'markersize',10)
plot(xAUCPart1,yAUCPart1,'k--','linewidth',1.5)
plot(xAUCPart2,yAUCPart2,'k--','linewidth',1.5)

set(gca,'fontsize',16)
xlabel('Protrusion Depth')
ylabel('Strength of Synapse (a.u.)')
xticks([-0.8:.2:0.8])
xticklabels({'-0.8','','-0.4','','0','','0.4','','0.8'})
ylim([0 18])
yticks([0:2:16])
yticklabels({'0','','4','','8','','14','','16'})

subplot(1,3,3); hold on;
plot((1-domain_x_length_vec(:)),half_max_width(:),'ko','linewidth',1.5,'markersize',10)
plot(xHMWPart1,yHMWPart1,'k--','linewidth',1.5)
plot(xHMWPart2,yHMWPart2,'k--','linewidth',1.5)

set(gca,'fontsize',16)
xlabel('Protrusion depth')
ylabel('Half-max width (a.u.)')
xticks([-0.8:.2:0.8])
xticklabels({'-0.8','','-0.4','','0','','0.4','','0.8'})
ylim([0 0.4])
yticks([0:0.05:0.4])
yticklabels({'0','','0.1','','0.2','','0.3','','0.4'})