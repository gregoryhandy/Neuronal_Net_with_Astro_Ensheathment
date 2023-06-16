%%
% This file contains all parameters, EXCEPT for prob_ensheathed and 
% ensh_strength, which are defined in the main function
%
% Key pamaters:
%   sigmaffwd: feedforward projection radius
%   sigmarec: recurrent projection radius
%
% Neuronal network result: sigmaffwd < sigmarec leads correlated balanced state
% Breaks down with astrocytes! Possible with sigmaffwd > sigmarec
%
%%
function [param]  = spatial_network_params_fn()

% Feedforward and recurrent
% connection widths
param.sigmaffwd= 0.10;
% For correlated state in the default neuronal network use sigmarec = 0.25
param.sigmarec = 0.05; %0.2;

% Number of neurons in network along each dimension 
% (e.g., the exc network is Ne1xNe1)
param.Ne1=200;
param.Ni1=100;
param.NF1=75;
param.Ne=param.Ne1*param.Ne1;
param.Ni=param.Ni1*param.Ni1;
param.NF=param.NF1*param.NF1;
param.N1=param.Ne1+param.Ni1;
param.N=param.Ne+param.Ni;

% Proportion of neurons that are excitatory
param.q=param.Ne/param.N;

% Outgoing connection probabilities scalings
param.pee0=0.05;
param.pei0=0.05;
param.pie0=0.05;
param.pii0=0.05;

% Connection strength scalings
param.jee=40;
param.jei=400;
param.jie=120; 
param.jii=400;

% average strength of connections (used in theory)
param.wee0=param.jee*param.pee0*param.q;
param.wei0=param.jei*param.pei0*(1-param.q);
param.wie0=param.jie*param.pie0*param.q;
param.wii0=param.jii*param.pii0*(1-param.q);

param.jeF=120;
param.jiF=120;
param.rF=.005;
param.peF0=.25;
param.piF0=.08;

param.weF0=param.jeF*param.peF0*(param.NF/param.N);
param.wiF0=param.jiF*param.piF0*(param.NF/param.N);

% Time constants
param.tausyne=6;
param.tausyni=5;
param.tausynF=param.tausyne;

% Random Excitatory neurons to record 
% currents and voltages from
param.nrecord0=200;
param.Irecord=randi(param.Ne1,2,param.nrecord0);

% Neuron params
param.gl=[1/15 1/10];
param.Cm=[1 1];
param.vlb=[-100 -100];
param.vth=[-10 -10];
param.DeltaT=[2 .5]; 
param.vT=[-50 -50]; %mV
param.vl=[-60 -60]; %mV
param.vre=[-65 -65];
param.tref=[1.5 .5];

% Width of recurrent connections
param.sigmaee=param.sigmarec;
param.sigmaie=param.sigmarec;
param.sigmaei=param.sigmarec;
param.sigmaii=param.sigmarec;
param.sigmaeF=param.sigmaffwd;
param.sigmaiF=param.sigmaffwd;

% Widths of connections in terms of number
% of postsynaptic neurons
param.betaee=param.sigmaee*(param.Ne1);
param.betaei=param.sigmaei*(param.Ne1);
param.betaie=param.sigmaie*(param.Ni1);
param.betaii=param.sigmaii*(param.Ni1);
param.betaeF=param.sigmaeF*(param.Ne1);
param.betaiF=param.sigmaiF*(param.Ni1);

% Number of outgoing connections
param.Kee=param.pee0*param.Ne;
param.Kei=param.pei0*param.Ne;

param.Kie=param.pie0*param.Ni;
param.Kii=param.pii0*param.Ni;

param.KeF=param.peF0*param.Ne;
param.KiF=param.piF0*param.Ni;

param.total_synapses = param.Ne*(param.Kee+param.Kie)+param.Ni*(param.Kei+param.Kii);

% Actual connection weights
param.Jee=(param.jee/sqrt(param.N));
param.Jei=-(param.jei/sqrt(param.N));

param.Jie=(param.jie/sqrt(param.N));
param.Jii=-(param.jii/sqrt(param.N));

param.JeF=param.jeF/sqrt(param.N);
param.JiF=param.jiF/sqrt(param.N);

%% Time stuff
% Length, time bin, and burn-in period for simulation
% T=22000;
param.T=5000;
param.dt=.1;
param.Nt=round(param.T/param.dt);
param.Tburn=2000;
param.nburn=round(param.Tburn/param.dt);

% Maximum number of spikes.
% Simulation will terminate with a warning if this is exceeded
param.maxns=param.N*param.T*.1;

end