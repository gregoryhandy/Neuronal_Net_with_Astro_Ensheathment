%%
% Plots all main figures for "Investigating the ability of astrocytes to 
% drive neural network synchrony"
% https://www.biorxiv.org/content/10.1101/2022.09.26.508928v1
%
% See the individual files for additional details 
%%

clear; close all; clc;

% DiRT simulations
plotFig1;

% Non-spatial, one population
plotFig3and4;

% Non-spatial, two populations
plotFig5;

% Spatial 
plotFig6; 

% Non-spatial, two population, both exc. and inh ensheathment
plotFig7;

% Non-spatial, one population supplemental figure
plotS1Fig;

% Spatial supplemental figure
plotS2Fig;