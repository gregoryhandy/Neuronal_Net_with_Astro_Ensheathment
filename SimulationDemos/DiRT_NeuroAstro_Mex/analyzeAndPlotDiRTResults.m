%%
% Analyze and plot the official results of the DiRT simulations (reproduces
% Fig 1)
% 
% Loading DiRT_Results adds three items to the workspace:
%   1) DiRT_results: # time points x 2 x # trials x # phi values
%       The first column contains the time vectors, while the second
%       contains the number of activate receptors
%   2) params: parameters used across all simulations
%       Except for domain_x_total and cleft_x_total which depend on the phi
%       values
%   3) phi_vec: vector of the phi values considered
%%
clear; close all; clc;

%% Load the data
load('./SimResults/DiRT_Results_official')

%% Ave the results and get curve statistics

max_val = zeros(length(phi_vec),1);
half_max_start= zeros(length(phi_vec),1);
half_max_end= zeros(length(phi_vec),1);
half_max_width= zeros(length(phi_vec),1);
AUC = zeros(length(phi_vec),1);
legend_names = cell(length(phi_vec),1);

for ii = 1:length(phi_vec)
    
    % Average across trials
    DiRT_results_ave = squeeze(mean(DiRT_results,3));
    
    % Half-max width
    [max_val(ii), index_start ] = max(DiRT_results_ave(:,2,ii));
    half_max_start(ii) = find(DiRT_results_ave(:,2,ii) > max_val(ii)/2,1);
    half_max_end(ii) = find(DiRT_results_ave(:,2,ii) < max_val(ii)/2 & ...
        DiRT_results_ave(:,1,ii)>DiRT_results_ave(half_max_start(ii),1,ii),1);
    half_max_width(ii) = DiRT_results_ave(half_max_end(ii),1,ii)-DiRT_results_ave(half_max_start(ii),1,ii);

    % Area under the curve (i.e., synaptic strength)
    AUC(ii) = trapz(DiRT_results_ave(:,1,ii),DiRT_results_ave(:,2,ii));

    legend_names{ii} = sprintf('%.1f',phi_vec(ii));  
end

%% Plot example activation time courses
color_scheme = [0.85 0.325 0.098; 0 0.447 0.741; 0.466 0.674 0.188; 0 0 0];
ex_indices = [9,13];

f = figure(113); clf; 
f.Position(3:4) = [1271 420];
subplot(1,3,1); hold on; h=[];
for ii = 1:length(ex_indices)
    
    h(ii) = plot(DiRT_results_ave(:,1,ex_indices(ii)),...
        DiRT_results_ave(:,2,ex_indices(ii)),'linewidth',3,'color',color_scheme(ii,:));
    plot(DiRT_results_ave(half_max_start(ex_indices(ii)),1,ex_indices(ii)), ...
        max_val(ex_indices(ii))/2,'.','markersize',20,'color',color_scheme(ii,:))
    plot(DiRT_results_ave(half_max_end(ex_indices(ii)),1,ex_indices(ii)), ...
        max_val(ex_indices(ii))/2,'.','markersize',20,'color',color_scheme(ii,:))
     
    t_width = [DiRT_results_ave(half_max_start(ex_indices(ii)),1,ex_indices(ii))...
        :0.01:DiRT_results_ave(half_max_end(ex_indices(ii)),1,ex_indices(ii))];
    plot(t_width,max_val(ex_indices(ii))/2*ones(length(t_width),1),'--',...
        'linewidth',2,'color',color_scheme(ii,:))
end
set(gca,'fontsize',16)
xlim([0 1])
legend(h,strcat('\phi= ',legend_names(ex_indices)))
xlabel('Time (a.u.)')
ylabel('# of Activated Receptors')
yticks([0:10:50])

%% Fit the piecewise linear model for the full range of parameter values

% Fit the model in the two regimes
% Saturated regime
xAUCPart1 = [phi_vec(1:7) -0.15];
yAUCPart1 = mean(AUC(1:6))*ones(8,1);

xHMWPart1 = [phi_vec(1:7) -0.15];
yHMWPart1 = mean(half_max_width(1:6))*ones(8,1);

% Fit the linear regime for the AUC
[linear_model_AUC,goodness_AUC] = fit(phi_vec(7:end)',AUC(7:end),'poly1');
xAUCPart2 = [-0.15:0.001:1];
fitStart = [xAUCPart1(end), yAUCPart1(end)];
yAUCPart2 = linear_model_AUC.p1*(xAUCPart2-fitStart(1))+fitStart(2);

% Fit the linear regime for the half-max width
[linear_model_HMW,~] = fit(phi_vec(7:end-4)',...
    half_max_width(7:end-4),'poly1');
xHMWPart2 = [-0.15:0.001:1];
fitStart = [xHMWPart1(end), yHMWPart1(end)];
yHMWPart2 = linear_model_HMW.p1*(xHMWPart2-fitStart(1))+fitStart(2);

%% Panels with data and piecewise linear fits
subplot(1,3,2); hold on;
plot(phi_vec,AUC,'ko','linewidth',1.5,'markersize',10)
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
plot(phi_vec,half_max_width,'ko','linewidth',1.5,'markersize',10)
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


