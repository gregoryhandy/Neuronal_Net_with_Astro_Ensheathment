%%
% Simulation results for the non-spatial one population model for networks
% with and without astrocyte ensheathment
%
% Creates S2 Fig
%%

% Load the data
load('./CompressedData/NonspatialOnePop/S2Fig_both_ex.mat')
load('./CompressedData/NonspatialOnePop/S2Fig_J_ex.mat')
load('./CompressedData/NonspatialOnePop/S2Fig_tau_ex.mat')

%% Create all panels
f=figure(9); clf; hold on
f.Position(3:4)=[1395 303];

% Create the raster plots
for ii = 1:3
    if ii == 1
        currData = S2Fig_both_ex;
    elseif ii == 2
        currData = S2Fig_J_ex;
    else
        currData = S2Fig_tau_ex;
    end
    
    subplot(1,4,ii); hold on;
    excNeurons = find(currData.sRaster(2,:)<100);
    inhNeurons = find(currData.sRaster(2,:)>=5000 & currData.sRaster(2,:)<=5100);
    
    plot(currData.sRaster(1,excNeurons)/1000,currData.sRaster(2,excNeurons),'.')
    plot(currData.sRaster(1,inhNeurons)/1000,currData.sRaster(2,inhNeurons)-4901,'.')
    xlim([0 2.5])
    set(gca,'fontsize',16)
    xlabel('Time (s)')
    ylabel('Neuron Index')
end

% Create the bar plot
subplot(1,4,4); 
bar([0 1],[S2Fig_J_ex.numSuccesses,S2Fig_tau_ex.numSuccesses]/10*100)
xticklabels({'J','tau'})
xlabel('Parameter')
ylabel('Percent Synchronous')
set(gca,'fontsize',16)
box off



