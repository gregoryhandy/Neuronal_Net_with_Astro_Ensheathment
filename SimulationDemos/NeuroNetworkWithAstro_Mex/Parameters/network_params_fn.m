%%
% This file contains all parameters, EXCEPT for num_ffwd_inputs,
% prob_ensheathed and ensh_strength, which are defined in the main 
% function
%%
function [param]  = network_params_fn()

% Number of neurons in each population
param.Ne=10000;
param.Ni=10000;
param.N=param.Ne+param.Ni;
param.Ne1=param.Ne/2;   % the subpopulation recent different background input
param.Ni1=param.Ni/2;
param.q=param.Ne/param.N;

% Outgoing connection probabilities 
param.pee0=.25;
param.pei0=.25;
param.pie0=.25;
param.pii0=.25;

% NOTE: these must be whole numbers, choose Ne, Ni appropriately 
% Number of outgoing connections
% So each excitatory neuro has exactly Kee + Kie outgoing connections
param.Kee=param.pee0*param.Ne;
param.Kei=param.pei0*param.Ne;
param.Kie=param.pie0*param.Ni;
param.Kii=param.pii0*param.Ni;

param.total_synapses = param.Ne*(param.Kee+param.Kie)+param.Ni*(param.Kei+param.Kii);

% Connection weight scalings
param.jee=12.5;
param.jei=50;
param.jie=20;
param.jii=50;

% Actual connection weights
param.Jee=(param.jee/sqrt(param.N));
param.Jei=-(param.jei/sqrt(param.N));
param.Jie=(param.jie/sqrt(param.N));
param.Jii=-(param.jii/sqrt(param.N));

% average strength of connections (used in theory)
param.wee0=param.jee*param.pee0*param.q;
param.wei0=param.jei*param.pei0*(1-param.q);
param.wie0=param.jie*param.pie0*param.q;
param.wii0=param.jii*param.pii0*(1-param.q);

% Determines synaptic decay and height (normalized to 1)
param.tausyne = 5;
param.tausyni = 4;

% Timescale, mean and std of ffwd input current
param.taux=40;
param.mxe=0.015*sqrt(param.N);
param.mxi=0.01*sqrt(param.N);
param.sigma_s = 0.1; 

param.mxe0=param.mxe/sqrt(param.N);
param.mxi0=param.mxi/sqrt(param.N);

% Pick the neuron indices for which to record currents
if(param.N>500)
    param.Nrecord=500;
else
    param.Nrecord = 5;
end

% using this to prevent duplications with the small network
% but this will be much slower for larger networks
% Irecord=sort(randperm(Ne,Nrecord)); 
param.Irecord=sort(randi(param.Ne,param.Nrecord,1));

%% Neuron params (units are in mV and ms)
param.gl=[1/15 1/10];    % leak conductance
param.C=[1 1];           % membrane capacitance
param.Vlb=[-100 -100];   % lower bound for cell potential
param.Vth=[-10 -10];     % spike threshold 
param.DeltaT=[2 0.5];    % tuning param for the?spike-generating current
param.VT=[-50 -50];      % soft spike threshold
param.Vleak=[-60 -60];   % leak potential
param.Vre=[-65 -65];     % reset threshold
param.tref=[1.5 .5];     % refractory time 

%% Time stuff
param.T=5000;               % total time of simulation
param.dt=0.05;                % time step
param.Nt=round(param.T/param.dt);        % number of time steps
param.Tburn=2000;            % timed allowed simulation to stabilize
param.nburn=round(param.Tburn/param.dt); 

% convert the refractory time to time steps
param.Ntref(1)=round(param.tref(1)/param.dt);
param.Ntref(2)=round(param.tref(2)/param.dt);

% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
param.maxns=param.N*param.T*.02;

end

