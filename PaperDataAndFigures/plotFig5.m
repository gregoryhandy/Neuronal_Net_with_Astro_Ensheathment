%%
% Simulation results for the non-spatial two populations model for networks
% with and without astrocyte ensheathment of exc synapses
%%

%% Figure 5A and B: correlation across neurons for a specific example

% Load the data
prob_ensheathed_I = 0.0;
prob_ensheathed_E = 0.8;
ensh_strength = 0.5;
num_ffwd_inputs=2;
num_trials=5;

filename_temp = ...
    sprintf('./CompressedData/NonspatialTwoPop/corr_data_%.2fprobEnI_%.2fprobEnE_%.2fenS_%dinputs_%dtrials.mat',...
    prob_ensheathed_I,0.00,ensh_strength,num_ffwd_inputs,num_trials);
default_ex = load(filename_temp);

filename_temp = ...
    sprintf('./CompressedData/NonspatialTwoPop/corr_data_%.2fprobEnI_%.2fprobEnE_%.2fenS_%dinputs_%dtrials.mat',...
    prob_ensheathed_I,prob_ensheathed_E,ensh_strength,num_ffwd_inputs,num_trials);
astro_ex = load(filename_temp);

figure(5); clf;
subplot(2,3,1); hold on;
plot(mean(default_ex.b_all),mean(default_ex.h_all),'k','linewidth',2)
plot(mean(astro_ex.b_all),mean(astro_ex.h_all),'k--','linewidth',2)
ylim([0 1.5])
yticks([0:0.5:1.5])
set(gca,'fontsize',16)
xlabel('Spike Count Correlations')
ylim([0 1.4])
yticks([0:0.2:1.4])
title('All Neurons')
legend('Default Network','Network w/ Astro')

subplot(2,3,2); hold on;
plot(mean(default_ex.b_opp),mean(default_ex.h_opp),'linewidth',2,'color',[0,128,128]/255)
plot(mean(default_ex.b_same),mean(default_ex.h_same),'linewidth',2,'color',[128,49,167]/255)
plot(mean(astro_ex.b_opp),mean(astro_ex.h_opp),'--','linewidth',2,'color',[0,128,128]/255)
plot(mean(astro_ex.b_same),mean(astro_ex.h_same),'--','linewidth',2,'color',[128,49,167]/255)
ylim([0 3])
yticks([0:1:3])
set(gca,'fontsize',16)
xlabel('Spike Count Correlations')
ylim([0 3])


%% Figure 5D: Create the raster plots for example simulations

% Load the data
load('./CompressedData/NonspatialTwoPop/astro_2Pop_ex.mat')
load('./CompressedData/NonspatialTwoPop/default_2Pop_ex.mat')

indicesPop1_d=find(default_2Pop_ex.sRaster(2,:)<=200);
indicesPop2_d=find(default_2Pop_ex.sRaster(2,:)>5000 & default_2Pop_ex.sRaster(2,:)<=5000+200);

indicesPop1_a=find(astro_2Pop_ex.sRaster(2,:)<=200);
indicesPop2_a=find(astro_2Pop_ex.sRaster(2,:)>5000 & astro_2Pop_ex.sRaster(2,:)<=5000+200);

subplot(2,3,4); hold on;
plot(default_2Pop_ex.sRaster(1,indicesPop1_d)/1000,default_2Pop_ex.sRaster(2,indicesPop1_d),...
    '.','markersize',10,'color',[171,35,40]/255)
plot(default_2Pop_ex.sRaster(1,indicesPop2_d)/1000,default_2Pop_ex.sRaster(2,indicesPop2_d)-4800,...
    '.','markersize',10,'color',[255,140,0]/255)
plot([2:0.1:3],200*ones(length([2:0.1:3]),1),'k-','linewidth',3)

title('Default Network')
xticks([2 2.5 3])
xlim([2 3])
yticks([100 300])
yticklabels({'Pop 1','Pop2'})
set(gca,'fontsize',16)
xlabel('Time (s)')

subplot(2,3,5); hold on;
plot(astro_2Pop_ex.sRaster(1,indicesPop1_a)/1000,astro_2Pop_ex.sRaster(2,indicesPop1_a),...
    '.','markersize',10,'color',[171,35,40]/255)
plot(astro_2Pop_ex.sRaster(1,indicesPop2_a)/1000,astro_2Pop_ex.sRaster(2,indicesPop2_a)-4800,...
    '.','markersize',10,'color',[255,140,0]/255)
plot([2:0.1:3],200*ones(length([2:0.1:3]),1),'k-','linewidth',3)
title('Network with Astrocytes')
xlim([2 3])
xticks([2 2.5 3])
yticks([100 300])
yticklabels({'Pop 1','Pop2'})
set(gca,'fontsize',16)
xlabel('Time (s)')

%% Figure 5C: correlation values across different levels of ensheathment

% Load the data
load('./CompressedData/NonspatialTwoPop/enEandIGridData_2Pop.mat')

% We only need the data for p_en^i = 0 for this plot
corrMatrix = squeeze(corrMatrix(1,:,:));
rateMatrixExc = squeeze(rateMatrixExc(1,:,:));

subplot(2,3,3); hold on;
surf(p_en,s_en,corrMatrix')
caxis([0.4 0.55])
colorbar
colormap(linspecer)
xlim([0.1 1])
ylim([0 0.5])
xlabel('Prob. of Ensheathment')
ylabel('Strength of Ensheathment')
title('Mean Correlations (Same Pop)')
set(gca,'fontsize',16)
box on

%% Figure 5E: excitatory firing rates vs. correlation panel
subplot(2,3,6); hold on;
colormap(linspecer)

for kk = 1:size(p_en,1)
    scatter(corrMatrix(:,kk),rateMatrixExc(:,kk),...
        200,p_en(kk,:).*s_en(kk,:)/100,'filled');
end
set(gca,'fontsize',16)
colorbar()

%Fit the data;
corrMatrixAdj = corrMatrix(:);
rateMatrixAdj = rateMatrixExc(:);

indices = isnan(corrMatrixAdj);
corrMatrixAdj(indices)=[];
rateMatrixAdj(indices)=[];

X_fit = [ones(length(corrMatrixAdj(:)),1) corrMatrixAdj(:)];
b_fit = X_fit\rateMatrixAdj(:);

x_temp = [min(corrMatrixAdj(:))*0.95:0.0001:max(corrMatrixAdj(:))*1.05];
plot(x_temp,b_fit(1) + b_fit(2)*x_temp,'k--','linewidth',1.5);
xlabel('Spike count correlations')
ylabel('Excitatory firing rate (Hz)')

model_vals = b_fit(1) + b_fit(2)*corrMatrixAdj(:);

r_squared = 1-sum((rateMatrixAdj(:)-model_vals).^2)/sum((rateMatrixAdj(:)-mean(rateMatrixAdj(:))).^2);
fprintf('Goodness of fit for Fig 5E: R^2 = %.2f \n',r_squared)
xlim([min(corrMatrixAdj(:))*0.95, max(corrMatrixAdj(:))*1.05])

% to save use:
% set(gcf,'renderer','painters')
% print('-dpdf','filename.pdf')

