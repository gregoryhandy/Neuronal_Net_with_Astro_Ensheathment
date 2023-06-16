%%
% Function Inputs
%   total_synapses: total number of synapses
%   distribution:  binomial with parameter p (prob. of success)
%       0 for not ensheathed
%       1 for ensheathed
%   d_parameter: distribution parameter
%       p (probability of selecting zero ensheathment
%   synapse_type (Note: IE reads, excitatory to inhibitory 
%       1: EE, 2: IE, 3: EI, and 4: II
%%
function [ensh_realizations,tausyn,J] =...
    create_ensheathment(N,total_synapses,prob_ensheathed_E,prob_ensheathed_I,...
    synapse_type,jee, jei, jie, jii, tausyne, tausyni, s_en)

    % flip a coin to determine if synapse is ensheathed  
    % Different coins!
    % Make sure the total number of synapses is even and the there are same
    % number of each type (can be generalized!)
    ensh_realizations = zeros(total_synapses,1);
    ensh_realizations(1:total_synapses/2) = binornd(1,prob_ensheathed_E,total_synapses/2,1);
    ensh_realizations((total_synapses/2+1):end) = binornd(1,prob_ensheathed_I,total_synapses/2,1);
    
    tausyn = zeros(4,1);
    J = zeros(total_synapses,1);
    
    % Actual connection weights (keep these fixed for now)
    Jee=(jee/sqrt(N)); Jei=-(jei/sqrt(N));
    Jie=(jie/sqrt(N)); Jii=-(jii/sqrt(N));
        
    % Synaptic kinetics take the form: 1/tau*exp(-t/tau)*Heaviside(t)
    % tausyn is a vector of taus for different synapses
    tausyn(1) = tausyne;               % default excitatory (tausyne0)
    tausyn(2) = tausyne*(1-s_en); % ensheathed excitatory (tausyne1)
    tausyn(3) = tausyni;               % default inhibitory (tausyni0)
    tausyn(4) = tausyni*(1-s_en); % ensheathed inhibitory (tausyni1)
        
    % Loop through the synapses and note the connection weight
    % Kept if in the future we want ensheahment to effect these
    for i = 1:total_synapses
        if synapse_type(i) == 1 
            J(i) = Jee;
        elseif synapse_type(i) == 2
            J(i) = Jie;
        elseif synapse_type(i) == 3
            J(i) = Jei;
        else
            J(i) = Jii;
        end
    end
    
end

