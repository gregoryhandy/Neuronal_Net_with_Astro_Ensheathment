%%
% Simulation results for the non-spatial two populations model for networks
% with and without astrocyte ensheathment of exc and inh synapses
%%

%% Figure 7A: Correlation surfaces across different values of p_en^i

% load the data
load('./CompressedData/NonspatialTwoPop/enEandIGridData_2Pop.mat')

% probability of ensheathment of inhibitory neurons
prob_ensheathed_I = 0:0.2:1.0;

figure(7); clf; 
for ii = 1:length(prob_ensheathed_I)
    if ii > 3
        subplot(2,4,ii+1); hold on;
    else
        subplot(2,4,ii); hold on;
    end
    surf(p_en,s_en,squeeze(corrMatrix(ii,:,:))')
    
    caxis([0.4 0.75])
    colorbar
    colormap(linspecer)
    xlim([0.1 1])
    title(sprintf('Prob. of Inh Ensheathment:%.2f',prob_ensheathed_I(ii)))
    box on
end

%% Figure 7B: Scatter plot relating firing rates with correlations
subplot(2,4,[4 8]); hold on;
colormap(linspecer)
for ii = 1:length(prob_ensheathed_I)
    for kk = 1:size(p_en,1)
        scatter(corrMatrix(ii,:,kk),rateMatrixExc(ii,:,kk),...
            200,p_en(kk,:).*s_en(kk,:)/100,'filled');
    end
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
fprintf('Goodness of fit for Fig 7B: R^2 = %.2f \n',r_squared)
xlim([min(corrMatrixAdj(:))*0.95, max(corrMatrixAdj(:))*1.05])

% to save use:
% set(gcf,'renderer','painters')
% print('-dpdf','filename.pdf')