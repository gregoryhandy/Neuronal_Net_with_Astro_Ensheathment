%%
% Simulation results for the non-spatial one population model for networks
% with and without astrocyte ensheathment
% Creates Figures 3 and 4
%%

%% Figure 3: Plot the raster and rolling average plots
% Load the data raster plot examples
load('./CompressedData/NonspatialOnePop/balancedCompressed.mat')
load('./CompressedData/NonspatialOnePop/brokenEx1Compressed.mat')
load('./CompressedData/NonspatialOnePop/brokenEx2Compressed.mat')

timeVec=balance.params.dt:balance.params.dt:balance.params.Nt;
tMin_plot = 1;
tvec_v2 = [0.05:0.05:5000]/1000;
numNeurons = 100;
excColor = [0.85 0.325 0.098];
inhColor = [0 0.447 0.741];

figure(3); clf;
for ii = 1:3
    if ii == 1
        currData = balance;
    elseif ii == 2
        currData = broken_ex1;
    elseif ii == 3
        currData = broken_ex2;
    end
    
    subplot(3,3,1+(ii-1)); hold on;
    plot(tvec_v2,currData.Ix1e,'k','linewidth',1.5)

    xlim([tMin_plot 4.5])
    set(gca,'fontsize',16)
    ylabel('Ffwd Input Current')
    box off
    
    subplot(3,3,4+(ii-1)); hold on;
    exc_neurons = find(currData.neuroninds <= numNeurons);
    inh_neurons = find(currData.neuroninds > numNeurons);
    plot(currData.s_raster(exc_neurons),currData.neuroninds(exc_neurons),...
        '.','color',excColor,'MarkerSize',5)
    plot(currData.s_raster(inh_neurons),currData.neuroninds(inh_neurons),...
        '.','color',inhColor,'MarkerSize',5)
    axis tight
    xlim([tMin_plot 4.5])
    set(gca,'fontsize',16)
    ylabel('Neuron index')
    box off
    
    subplot(3,3,7+(ii-1)); hold on;
    plot(currData.tvec/1000,currData.reSim_roll,'color',excColor,'linewidth',1.5)
    plot(currData.tvec/1000,currData.riSim_roll,'color',inhColor,'linewidth',1.5)
    xlim([tMin_plot 4.5])
    if ii == 2 || ii == 3 
        ylim([0 30])
    else
        ylim([2 7])
        yticklabels({'2','','4','','6'})
    end
    set(gca,'fontsize',16)
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    box off
end


%% Figure 4A: Plot the recurrent, feedforward, total input curves
Tmin=1000; Tmax=1850;
tIndices=find(timeVec>=Tmin & timeVec<=Tmax);

% Time window to average over for baseline
tAveIndices = find(timeVec>=Tmin & timeVec<=1800);

figure(4); clf; 
subplot(3,3,1); hold on; clearvars h; 
for ii = 1:2
    if ii == 1
        currData = balance; lineType = '-';
    else
        currData = broken_ex1; lineType='--';
    end
    
    baselineRec=mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices));
    baselineAll = mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices)+IF0_lowpass(tAveIndices));
    
    h(1+2*(ii-1)) = plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)...
        -baselineRec)/currData.Grheobase,lineType,'color',[0.55 0 0.55],'LineWidth',2);
    h(2+2*(ii-1)) = plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)+IF0_lowpass(tIndices)...
        -baselineAll)/currData.Grheobase,lineType,'color',[0.5 0.5 0.5],'LineWidth',2);
end

baselineFfwd = mean(IF0_lowpass(tAveIndices));
h(5) = plot(timeVec(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/currData.Grheobase,...
    'color','k','LineWidth',2);

legend([h(1),h(2),h(5)],{'Recurrent','Feedforward','Total'})
set(gca,'fontsize',16)
ylabel('Normalized Shared Current')
xlabel('Time (s)')
ylim([-2 2])
xlim([1000 1850]/1000)
xticks([1:0.2:1.9])
box off

%% Figure 4B: Plot the zoom-in panels
Tmin=500; Tmax=4000;
time=balance.params.dt:balance.params.dt:balance.params.Nt;
tIndices=find(time>=Tmin & time<=Tmax);

excColor = [0.85 0.325 0.098];
inhColor = [0 0.447 0.741];

% figure(3); clf; 
subplot(3,3,2); hold on;
plot(time(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/balance.Grheobase,'-','color','k','linewidth',1.5)

subplot(3,3,3); hold on;
plot(time(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/balance.Grheobase,'-','color','k','linewidth',1.5)

for ii = 1:2
    if ii == 1
        currData = balance; lineType = '-';
    else
        currData = broken_ex1; lineType='--';
    end
    
    baselineExc=mean(currData.IeAve(tAveIndices));
    baselineInh=mean(currData.IiAve(tAveIndices));
    baselineRec=mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices));
    baselineAll = mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices)+IF0_lowpass(tAveIndices));
    
    for jj = 1:2
        subplot(3,3,2+(jj-1)); hold on;
        plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)...
            -baselineRec)/currData.Grheobase,lineType,'color',[0.55 0 0.55],'LineWidth',2);
        plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)+IF0_lowpass(tIndices)...
            -baselineAll)/currData.Grheobase,lineType,'color',[0.5 0.5 0.5],'LineWidth',2);
        ylim([-2 2])
        
        subplot(3,3,5+(jj-1)); hold on;
        plot(time(tIndices)/1000,(currData.IeAve(tIndices)-baselineExc)/currData.Grheobase,lineType,...
            'color',excColor,'linewidth',1.5)
        ylim([-4 4])
        
        subplot(3,3,8+(jj-1)); hold on;
        plot(time(tIndices)/1000,(currData.IiAve(tIndices)-baselineInh)/currData.Grheobase,lineType,...
            'color',inhColor,'linewidth',1.5)
        ylim([-10 10])
    end
end

for ii = 1:3
    subplot(3,3,2+3*(ii-1))
    xlim([1 1.22])
    set(gca,'fontsize',16)
    
    subplot(3,3,3+3*(ii-1))
    xlim([1.66 1.88])
    set(gca,'fontsize',16)
end

subplot(3,3,8)
xlabel('Time (s)')

subplot(3,3,9)
xlabel('Time (s)')

legend('Default','Network w/ astrocytes')

%% Fig 4C: chance of synchrony for different levels of ensheathment

% Load the data
load('./CompressedData/NonspatialOnePop/syncGridData.mat')
 
subplot(3,3,4); hold on;
colormap(linspecer)
for i = 1:size(s_en,1)
    scatter(p_en(i,:),s_en(i,:),75,percentSynced(:,i)*100,'filled')
end
set(gca,'fontsize',16)
ylim([0.2 0.8])
xlim([0.2 1])
colorbar
xlabel('Prob. of Ensheathment')
ylabel('Strength of Ensheathment')

%% Fig 4D: mean field approx panel

mfApprox = load('./CompressedData/NonspatialOnePop/syncGridDataMeanFieldApprox');

subplot(3,3,7); hold on;
colormap(linspecer)
for i = 1:size(mfApprox.s_en,1)
    scatter(mfApprox.p_en(i,i),mfApprox.s_en(i,i),200,mfApprox.percentSynced(i,i),'filled')
end
set(gca,'fontsize',16)
yticks([.45 .55 .65 .75])
xlim([0.5 1])


percent_range_temp = [0.5:0.01:1];
plot(percent_range_temp,1-(3./(5*percent_range_temp)-(1-percent_range_temp)./percent_range_temp),'k--')
xlabel('Prob. of Ensheathment')
ylabel('Strength of Ensheathment')
colorbar
