%%
% Simulation results for the non-spatial one population model for networks
% with and without astrocyte ensheathment
%
% Creates S1 Fig for realization 2 from Fig 3C
%%

% Load the data
load('./CompressedData/NonspatialOnePop/balancedCompressed.mat')
load('./CompressedData/NonspatialOnePop/brokenEx2Compressed.mat')

figure(8); clf;

%% S1 Fig A: Plot the recurrent, feedforward, total input curves
Tmin=500; Tmax=4000;
timeVec=balance.params.dt:balance.params.dt:balance.params.Nt;
tIndices=find(timeVec>=Tmin & timeVec<=Tmax);

% Time window to average over for baseline
tAveIndices = find(timeVec>=Tmin & timeVec<=2700);

subplot(3,4,1); hold on; clearvars h; 
for ii = 1:2
    if ii == 1
        currData = balance; lineType = '-';
    else
        currData = broken_ex2; lineType='--';
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
xlim([1000 3000]/1000)
xticks([1:0.5:3])
box off

%% S1 Fig Panels B and C: Create the rastor plots
tMin_plot = 1;
tvec_v2 = [0.05:0.05:5000]/1000;
numNeurons = 100;
excColor = [0.85 0.325 0.098];
inhColor = [0 0.447 0.741];

for ii = 1:2
    if ii == 1
        currData = balance;
    elseif ii == 2
        currData = broken_ex2;
    end
    
    subplot(3,4,5+4*(ii-1)); hold on;
    exc_neurons = find(currData.neuroninds <= numNeurons);
    inh_neurons = find(currData.neuroninds > numNeurons);
    plot(currData.s_raster(exc_neurons),currData.neuroninds(exc_neurons),...
        '.','color',excColor,'MarkerSize',5)
    plot(currData.s_raster(inh_neurons),currData.neuroninds(inh_neurons),...
        '.','color',inhColor,'MarkerSize',5)
    axis tight
    xlim([tMin_plot 3.2])
    set(gca,'fontsize',16)
    ylabel('Neuron index')
    xlabel('Time (s)')
    box off 
end

%% S1 Fig D: Plot the zoom-in panels
Tmin=500; Tmax=4000;
time=balance.params.dt:balance.params.dt:balance.params.Nt;
tIndices=find(time>=Tmin & time<=Tmax);

excColor = [0.85 0.325 0.098];
inhColor = [0 0.447 0.741];


subplot(3,4,2); hold on;
plot(time(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/balance.Grheobase,'-','color','k','linewidth',1.5)

subplot(3,4,3); hold on;
plot(time(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/balance.Grheobase,'-','color','k','linewidth',1.5)

subplot(3,4,4); hold on;
plot(time(tIndices)/1000,(IF0_lowpass(tIndices)-baselineFfwd)/balance.Grheobase,'-','color','k','linewidth',1.5)

for ii = 1:2
    if ii == 1
        currData = balance; lineType = '-';
    else
        currData = broken_ex2; lineType='--';
    end
    
    baselineExc=mean(currData.IeAve(tAveIndices));
    baselineInh=mean(currData.IiAve(tAveIndices));
    baselineRec=mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices));
    baselineAll = mean(currData.Ie0_lowpass(tAveIndices)+currData.Ii0_lowpass(tAveIndices)+IF0_lowpass(tAveIndices));
    
    for jj = 1:3
        subplot(3,4,2+(jj-1)); hold on;
        plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)...
            -baselineRec)/balance.Grheobase,lineType,'color',[0.55 0 0.55],'LineWidth',2);
        plot(timeVec(tIndices)/1000,(currData.Ie0_lowpass(tIndices)+currData.Ii0_lowpass(tIndices)+IF0_lowpass(tIndices)...
            -baselineAll)/balance.Grheobase,lineType,'color',[0.5 0.5 0.5],'LineWidth',2);
        ylim([-2 2])
        
        subplot(3,4,6+(jj-1)); hold on;
        plot(time(tIndices)/1000,(currData.IeAve(tIndices)-baselineExc)/balance.Grheobase,lineType,...
            'color',excColor,'linewidth',1.5)
        ylim([-4 4])
        
        subplot(3,4,10+(jj-1)); hold on;
        plot(time(tIndices)/1000,(currData.IiAve(tIndices)-baselineInh)/balance.Grheobase,lineType,...
            'color',inhColor,'linewidth',1.5)
        ylim([-10 10])
    end
end

for ii = 1:3
    subplot(3,4,2+4*(ii-1))
    xlim([1 1.22])
    set(gca,'fontsize',16)
    
    subplot(3,4,3+4*(ii-1))
    xlim([1.66 1.88])
    set(gca,'fontsize',16)
    
    subplot(3,4,4+4*(ii-1))
    xlim([2.78 3])
    set(gca,'fontsize',16)
end

subplot(3,4,10)
xlabel('Time (s)')
subplot(3,4,11)
xlabel('Time (s)')
subplot(3,4,12)
xlabel('Time (s)')
subplot(3,4,2)
ylabel('Normalized Shared Current')

subplot(3,4,6)
ylabel('Normalized Exc. Current')

subplot(3,4,10)
ylabel('Normalized Inh. Current')

legend('Default','Network w/ astrocytes')