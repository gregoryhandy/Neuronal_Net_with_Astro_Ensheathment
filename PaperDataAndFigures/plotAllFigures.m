%%
% Plots all main figures for Handy and Borisyuk. Investigating the 
% ability of astrocytes to drive neural network synchrony. 2023
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

% Non-spatial, one population supplemental figure
% Ensheathment effects only tau or J
plotS2Fig;

% Spatial supplemental figure
plotS3Fig;