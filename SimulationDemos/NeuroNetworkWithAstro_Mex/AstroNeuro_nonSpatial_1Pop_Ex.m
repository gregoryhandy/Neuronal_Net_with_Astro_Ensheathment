%%
% This code runs multiple examples of the non-spatial network used in Handy
% and Borisyuk, 2023 with one feedforward input
%
% The examples show the default network (no astrocyte ensheathmenth) and
% three examples of the network with astrocyte ensheathment. The first two
% show synchronous behavior, and the last maintains asynchronous firing
% (similar to the examples shown in Fig. 3)
%
% Important notes:
%   The examples are run in parallel (via parfor)
%   The spiking code interfaces with mex c-code
%
% The simulation can be saved, but the .mat file will be rather large
% (>2 GBs; depends on length of trial)
% 
% Code takes  ~100 seconds (depending on the computer)
%%
clear; clc; close all;

% Add that folder plus all subfolders to the path.
restoredefaultpath;
folder = fileparts(which('AstroNeuro_nonSpatial_1Pop_Ex.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex file
mex ./Mex_Functions/AstroNeuro_nonSpatial_mex.c

%% Load parameters and set the number of pops and level of ensheathment

params = network_params_fn();

% Number of feedforward inputs
num_ffwd_inputs = 1;

% Note: a number in the range [0 1], where 0 means no ensheathment
s_en = 0.5;

%% Create feedforward input (fixed across all examples)
rng(923326); % seed the random number generator for consistency
[Ix1e, Ix1i, Ix2e, Ix2i] = ffwd_smooth_input(num_ffwd_inputs,...
    params.taux,params.dt,params.Nt,params.mxe,params.mxi,...
    params.sigma_s);

%% Create the connection matrix (fixed across all examples)
rng(232); % seed the random number generator for consistency
[target_neuron,firing_neuron,id_start_firing,...
    synapse_num_default,synapse_type] =...
    create_nonSpatial_network(params.Ne, params.Ni, params.Kee,...
    params.Kie, params.Kei, params.Kii);

%% Create vector of ensheathment values for each synapse
% This varies for each example 
% Each example is seeded for consistency
totalExamples = 4;
ensh_realizations = zeros(params.total_synapses,totalExamples);
fprintf('Choosing which synapses to ensheath\n')
for ii = 1:totalExamples
    % Default probabilities for ensheathment used for all examples
    prob_ensheathed_I = 0;
    prob_ensheathed_E = 0.7;
    if ii == 1 % default case (i.e., no ensheathment)
        rng(1032);
        prob_ensheathed_E = 0;
    elseif ii == 2
        rng(2124);
    elseif ii == 3
        rng(1033); 
    else
        rng(1723);
    end
    
    [ensh_realizations(:,ii),tausyn,J] = ...
        create_ensheathment(params.N, params.total_synapses,...
        prob_ensheathed_E,prob_ensheathed_I,...
        synapse_type, params.jee, params.jei, params.jie, params.jii,...
        params.tausyne, params.tausyni, s_en);
end

%% Clean up the indexing for the MEX code
% (i.e., offset the index by one for the MEX code (C starts indexing at 0))
Irecord_mex = int32(params.Irecord-1);
id_start_firing_mex = int32(id_start_firing-1);
firing_neuron_mex = int32(firing_neuron-1);
synapse_num_default_mex = int32(synapse_num_default-1);
target_neuron_mex = int32(target_neuron-1);

%% Random membrane potential initial conditions (same across all examples)
rng(131); % seed the random number generator for consistency
V0min=params.Vre(1);
V0max=params.VT(1);
V0=(V0max-V0min).*rand(params.N,1)+V0min;
v=V0;

%% Preallocate memory 
s = zeros(2, params.maxns,totalExamples);
Ie = zeros(params.Nrecord,params.Nt,totalExamples);
Ii = zeros(params.Nrecord,params.Nt,totalExamples);
vr = zeros(params.Nrecord,params.Nt,totalExamples);

%% Run the mex code
tic;
fprintf('Running the spiking simulations (may take a few minutes)\n')
parfor ii=1:totalExamples
    [s(:,:,ii),Ie(:,:,ii),Ii(:,:,ii),vr(:,:,ii)]=AstroNeuro_nonSpatial_mex(Ix1e,Ix2e,Ix1i,Ix2i,...
        params.Ne,params.Ni,params.Ne1,params.Ni1,...
        J,params.C,params.gl,params.Vleak,params.DeltaT,params.VT,...
        params.tref,params.Vth,params.Vre,params.Vlb,...
        tausyn,V0,params.T,params.dt,params.maxns,Irecord_mex,...
        id_start_firing_mex, firing_neuron_mex, synapse_num_default_mex,...
        target_neuron_mex,ensh_realizations(:,ii),s_en);    
end
toc;

%% Get the rolling averages for the firing rates
sliding_window = 100;
tvec = [0:sliding_window:5000];
reSim_roll = zeros(length(tvec),totalExamples);
riSim_roll = zeros(length(tvec),totalExamples);
for ii = 1:totalExamples
    for tt = 1:length(tvec)
        reSim_roll(tt,ii) = 1000*nnz(s(1,:,ii)>sliding_window*(tt-1) & ...
            s(1,:,ii)<sliding_window*tt & s(2,:,ii)<=params.Ne)/(params.Ne*(sliding_window));
        riSim_roll(tt,ii) = 1000*nnz(s(1,:,ii)>sliding_window*(tt-1) & ...
            s(1,:,ii)<sliding_window*tt & s(2,:,ii)>params.Ne)/(params.Ni*(sliding_window));
    end
end

%% Plot the results
colorScheme = [0.85 0.325 0.098;0 0.447 0.741];
Tmin = 1000; Tmax = 4500;
mrkrsz=5; % Size of circles plotted
nplot=100; % Number of cells to plot from each population

f=figure(1); clf;
f.Position(3:4) =[1133 693];
for ii = 1:totalExamples
    
    % Plot the spike times (i.e., raster plots)
    subplot(2,4,ii); hold on;
    % indices to plot from exc population
    Iplot1=find(s(1,:,ii)>=Tmin & s(1,:,ii)<=Tmax & s(2,:,ii)<=nplot);
    % indices to plot from inh population
    Iplot2=find(s(1,:,ii)>=Tmin & s(1,:,ii)<=Tmax & s(2,:,ii)>params.Ne & s(2,:,ii)<=params.Ne+nplot);
    plot((s(1,Iplot1,ii))/1000,s(2,Iplot1,ii),'.','color',colorScheme(1,:),'MarkerSize',mrkrsz)
    plot((s(1,Iplot2,ii))/1000,s(2,Iplot2,ii)-params.Ne+nplot,'.','color',colorScheme(2,:),'MarkerSize',mrkrsz)
    xticks([1:4])
    xlim([Tmin Tmax]/1000)
    set(gca,'fontsize',16)
    ylabel('Neuron Index') 
    xlabel('Time (s)')

    % Plot the rolling spike times
    subplot(2,4,ii+4); hold on;
    plot(tvec/1000,reSim_roll(:,ii),'color',colorScheme(1,:),'linewidth',1.5)
    plot(tvec/1000,riSim_roll(:,ii),'color',colorScheme(2,:),'linewidth',1.5)
    xlim([Tmin Tmax]/1000)
    xticks([1:4])
    xlim([Tmin Tmax]/1000)
    set(gca,'fontsize',16)
    ylabel('Firing rate (Hz)')
    xlabel('Time (s)')
    if ii == 1
        legend('Exc.','Inh.')
        ylim([2 10])
    elseif ii ==4 
        ylim([2 10])
    end
end

for ii = 1:totalExamples
    subplot(2,4,ii);
    if ii == 1
        title('Default network')
    else
        title(sprintf('Network w/ astro \n(Ex. %d)',ii-1));
    end
end
