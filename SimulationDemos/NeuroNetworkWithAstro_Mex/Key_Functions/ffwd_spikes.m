%%
% Generate spike times for the ffwd network
%%
function [ sF ] = ffwd_spikes(T,rF,NF,NF1)


% rF*NF*T is the mean number of spikes that will occur
tempspikes=sort(T*rand(poissrnd(rF*NF*T),1));

% Store in the data format expected by the C code,
% where neuron indices are assigned randomly
sF=zeros(3,numel(tempspikes));    % preallocate, numel: number of elements
sF(1,:)=tempspikes;
clear tempspikes;
sF(2,:)=ceil(rand(1,size(sF,2))*NF1);  % x coordinate of firing neuron
sF(3,:)=ceil(rand(1,size(sF,2))*NF1);  % y coordinate of firing neuron

end

