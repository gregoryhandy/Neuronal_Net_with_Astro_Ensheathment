%%
% Simulation results for the spatial model for networks with and without 
% astrocyte ensheathment of exc synapses
%%

% Default network with no spatial correlations
load('./CompressedData/Spatial/default_net_ex1'); load('./CompressedData/Spatial/astro_net_ex1');  

% Default network with spatial correlations
load('./CompressedData/Spatial/default_net_ex2'); load('./CompressedData/Spatial/astro_net_ex2');

%% Figure 6A: Spatial raster plots
% Note: the networks without default spatial correlations

figure(6); clf;

% Starting time for the time windows (in msec)
% Time '0' is mapped to the first number in this vector
timeStamps = 1000+ [64 84 104];

% Loop across the datasets
for jj = 1:2
    if jj == 1
        currData = default_net_ex1;
    else
        currData = astro_net_ex1;
    end
    
    %Loop across the time stamps
    for ii = 1:3
        
        % We consider a 3ms time window for the spikes
        Iplot=find(currData.s(1,:)>=timeStamps(ii) &...
            currData.s(1,:)<=(timeStamps(ii)+3));
        
        neuroninds_1=currData.s(2,Iplot);
        neuroninds_2=currData.s(3,Iplot);
        
        subplot(4,2,jj+2*(ii-1))
        plot(neuroninds_1,neuroninds_2,'k.','MarkerSize',10)
        axis([0 200 0 200])
        temp_title = sprintf('Time %d ms',20*(ii-1));
        if ii == 1 && jj == 1
            title(strcat(temp_title, ' (default)'))
        elseif ii == 1 && jj == 2
            title(strcat(temp_title, ' (astro)'))
        else
            title(temp_title);
        end
        set(gca,'fontsize',16)
        xticklabels([])
        yticklabels([])
    end
end

%% Figure 6B: Correlation curves across space

color_scheme = [0 0 0; 1 0.5 0];
for ii = 1:2
    
    if ii == 1
        currDataAstro = astro_net_ex1;
        currDataDeafult = default_net_ex1;
    else
        currDataAstro = astro_net_ex2;
        currDataDeafult = default_net_ex2;
    end
    
    
    subplot(4,2,6+ii); hold on;
    h2 = plot(currDataDeafult.distbins,currDataDeafult.SCcorrSimMean,...
        'color',color_scheme(1,:),'Linewidth',2);
    h1 = plot(currDataAstro.distbins,currDataAstro.SCcorrSimMean,...
        'color',color_scheme(2,:),'Linewidth',2);
    
    plot(0:.1:1,zeros(length(0:.1:1),1),'k--','linewidth',1)
    set(gca,'fontsize',16)
    xlabel('Distance (a.u.)')
    ylabel('Spike count correlations')
    xticks([0 0.25 0.5])
    
    xlim([-.02 .7])
    box off
    
    if ii == 1
        title('\alpha_{rec} < \alpha_{ffwd}')
        legend([h2 h1], {'Default','Network w/ Astro'})
    else
        title('\alpha_{rec} > \alpha_{ffwd}')
    end
end
