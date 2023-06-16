%%
% This code runs multiple examples of the non-spatial network used in Handy
% and Borisyuk, 2023 with two feedforward input
%
% The examples show the default network (no astrocyte ensheathmenth) and
% one examples of the network with astrocyte ensheathment. Similar to Fig 5
% in the main text
%
% Important notes:
%   The examples are run in parallel (via parfor)
%   The spiking code interfaces with mex c-code
%
% The simulation can be saved, but the .mat file will be rather large
% (>2 GBs; depends on length of trial)
% 
% Code takes  ~60 seconds (depending on the computer)
%%
clear; clc; close all;

% Add that folder plus all subfolders to the path.
restoredefaultpath;
folder = fileparts(which('AstroNeuro_nonSpatial_2Pop_Ex.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex file
mex ./Mex_Functions/AstroNeuro_nonSpatial_mex.c

%% Load parameters and set the number of pops and level of ensheathment
params = network_params_fn();

% Number of feedforward inputs
% 1 for async default, 2 for sync default
num_ffwd_inputs = 2;

% probability and strength of ensheathment
prob_ensheathed_I = 0;
prob_ensheathed_E = 0.8;

% Note: a number in the range [0 1], where 0 means no ensheathment
s_en = 0.5;

%% Create connection matrix, initial conditions, etc.

% Fix the ffwd inputs and connectivity across trials
randSeed = 123;
rng(randSeed,'twister');

% Feedforward inputs
[Ix1e, Ix1i, Ix2e, Ix2i] = ffwd_smooth_input(num_ffwd_inputs,...
    params.taux,params.dt,params.Nt,params.mxe,params.mxi,...
    params.sigma_s);

% Create the connection matrix
[target_neuron,firing_neuron,id_start_firing,...
    synapse_num_default,synapse_type] =...
    create_nonSpatial_network(params.Ne, params.Ni, params.Kee,...
    params.Kie, params.Kei, params.Kii);

% Create vector of ensheathment values for each synapse
ensh_realizations = zeros(params.total_synapses,2);

% Note: only example two has astrocyte ensheathment
[ensh_realizations(:,2),tausyn,J] = ...
    create_ensheathment(params.N, params.total_synapses,...
    prob_ensheathed_E,prob_ensheathed_I,...
    synapse_type, params.jee, params.jei, params.jie, params.jii,...
    params.tausyne, params.tausyni, s_en);

% offset the index by one for the MEX code (C starts indexing at 0)
Irecord_mex = int32(params.Irecord-1);
id_start_firing_mex = int32(id_start_firing-1);
firing_neuron_mex = int32(firing_neuron-1);
synapse_num_default_mex = int32(synapse_num_default-1);
target_neuron_mex = int32(target_neuron-1);

% Random membrane potential initial conditions
V0min=params.Vre(1);
V0max=params.VT(1);
V0=(V0max-V0min).*rand(params.N,1)+V0min;
v=V0;

%% Preallocate memory 
totalExamples = 2;
s = zeros(2, params.maxns,totalExamples);
Ie = zeros(params.Nrecord,params.Nt,totalExamples);
Ii = zeros(params.Nrecord,params.Nt,totalExamples);
vr = zeros(params.Nrecord,params.Nt,totalExamples);

%% Run the mex code
tic;
fprintf('Running the simulation\n')
parfor ii = 1:totalExamples
[s(:,:,ii),Ie(:,:,ii),Ii(:,:,ii),vr(:,:,ii)]=AstroNeuro_nonSpatial_mex(Ix1e,Ix2e,Ix1i,Ix2i,...
    params.Ne,params.Ni,params.Ne1,params.Ni1,...
    J,params.C,params.gl,params.Vleak,params.DeltaT,params.VT,...
    params.tref,params.Vth,params.Vre,params.Vlb,...
    tausyn,V0,params.T,params.dt,params.maxns,Irecord_mex,...
    id_start_firing_mex, firing_neuron_mex, synapse_num_default_mex,...
    target_neuron_mex,ensh_realizations(:,ii),s_en);
end
toc;

%% Plot the results
colorScheme = [0.85 0.325 0.098;1, 0.549, 0];
Tmin = 3000; Tmax = 4000;
mrkrsz=10; % Size of circles plotted
nplot = 100;
f=figure(1); clf;
f.Position(3:4) =[505 297];
for ii = 1:totalExamples
    
    % Plot the spike times (i.e., raster plots)
    subplot(1,2,ii); hold on;
    % indices to plot from exc population
    Iplot1=find(s(1,:,ii)>=Tmin & s(1,:,ii)<=Tmax & s(2,:,ii)<=nplot);
    % indices to plot from other exc population
    Iplot2=find(s(1,:,ii)>=Tmin & s(1,:,ii)<=Tmax & s(2,:,ii)>params.Ne/2 & s(2,:,ii)<=params.Ne/2+nplot);
    plot((s(1,Iplot1,ii))/1000,s(2,Iplot1,ii),'.','color',colorScheme(1,:),'MarkerSize',mrkrsz)
    plot((s(1,Iplot2,ii))/1000,s(2,Iplot2,ii)-params.Ne/2+nplot,'.','color',colorScheme(2,:),'MarkerSize',mrkrsz)
    plot([Tmin Tmax]/1000, nplot*ones(2,1),'k','linewidth',1.5)
    xticks([1:.5:4])
    xlim([Tmin Tmax]/1000)
    set(gca,'fontsize',16)
    ylabel('Neuron Index') 
    xlabel('Time (s)')
    
    if ii == 1
        title('Default Network')
    else
        title('Network w/ Astro')
    end
end
