%%
% This code runs multiple examples of the spatial network used in Handy
% and Borisyuk, 2023
%
% The examples show the default network (no astrocyte ensheathmenth) and
% one example of the network with astrocyte ensheathment. Similar to Fig 6
% in the main text
%
% Important notes:
%   The examples are run in parallel (via parfor)
%   The spiking code interfaces with mex c-code
%
% The simulation can be saved, but the .mat file will be rather large
% (>2 GBs; depends on length of trial)
% 
% Code takes  ~4 minutes (depending on the computer) for 5 seconds of
% simulation time
%   Longer than the non-spatial code since there are significantly more
%   neurons
%   Should run for much longer (in simulation time) if estimating the 
%   spatial correlations
%
%%
clear; close all; clc;

% Add all subfolders to path
% Determine where your m-file's folder is.
restoredefaultpath;
folder = fileparts(which('AstroNeuro_Spatial_Ex.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex file
mex ./Mex_Functions/AstroNeuro_Spatial_mex.c

%% Load the parameters
[params]  = spatial_network_params_fn();

% probability of ensheathment
prob_ensheathed_E = 0.9;
prob_ensheathed_I = 0;
% ensheathment strength
% 1 means no difference, 0 kills the response
s_en = 0.5;

%% Create connection matrix, initial conditions, etc.
rng(790); % seed the random number generator for consistency

% feedfoward input spikes
sF = ffwd_spikes(params.T, params.rF, params.NF, params.NF1);

% create the connection matrix
% Note: indices already offset by 1 for MEX
fprintf('Creating the connection matrix\n')
[Wex1, Wex2, Wix1, Wix2,target_neuron,firing_neuron,...
    id_start_firing,synapse_num_default,synapse_type] = ...
    create_spatial_network(params.Ne1, params.Kee, params.betaee, params.Kie, params.betaie,...
    params.Ni1, params.Kii, params.betaii, params.Kei, params.betaei,...
    params.NF1, params.KeF, params.betaeF, params.KiF, params.betaiF);

fprintf('Choosing which synapses to ensheath\n')
% Create vector of ensheathment values for each synapse 
ensh_realizations = zeros(params.total_synapses,2);
[ensh_realizations(:,2),tausyn,J] = create_ensheathment(params.N, ...
    params.total_synapses, prob_ensheathed_E,prob_ensheathed_I,...
    synapse_type, params.jee, params.jei, params.jie, params.jii,...
    params.tausyne, params.tausyni, s_en);


% Random initial membrane potentials
V0min=params.vre(1);
V0max=params.vT(1);
V0=(V0max-V0min).*rand(params.N,1)+V0min;

%% Preallocate memory
s = zeros(3,params.maxns,2);

%% Run the mex code (make sure it is already compiled)
fprintf('Running the simulations\n')
tic;
parfor ii = 1:2
[s(:,:,ii),~,~,~,~]=AstroNeuro_Spatial_mex(sF,params.NF1,params.Ne1,params.Ni1,params.JeF,params.JiF,... % 1-6 inputs
    J,params.Cm,params.gl,params.vl,params.DeltaT,params.vT,params.tref,params.vth,params.vre,params.vlb,...                % 7-16
    tausyn,s_en,params.tausynF,V0,params.T,params.dt,params.maxns,params.Irecord,...     % 17-24
    id_start_firing, firing_neuron, synapse_num_default,...  % 25-27
    target_neuron,ensh_realizations(:,ii), ... % 28-29
    params.KeF,params.KiF,Wex1, Wex2, Wix1, Wix2);     % 30-35
end
toc;

%% Create the raster plots
Tmin=1000;
f = figure(1); clf;
f.Position([3 4]) = [858 528];
for ii = 1:2
    for tt = 1:3
        subplot(2,3,tt+3*(ii-1))
        Iplot=find(s(1,:,ii)>=(Tmin+20*(tt-1)) & s(1,:,ii)<=(Tmin+20*(tt-1)+4));
        
        neuroninds_1=s(2,Iplot,ii);
        neuroninds_2=s(3,Iplot,ii);
        
        figure(1)
        plot(neuroninds_1,neuroninds_2,'k.','MarkerSize',10)
        axis([0 200 0 200])
        current_time = Tmin+20*(tt-1);
        temp_title = sprintf('Time %d ms',current_time-Tmin);
        title(temp_title)
        set(gca,'fontsize',16)
        xticklabels([])
        yticklabels([])
    end
end
subplot(2,3,1)
ylabel('Default network')
subplot(2,3,4)
ylabel('Network w/ astro')

