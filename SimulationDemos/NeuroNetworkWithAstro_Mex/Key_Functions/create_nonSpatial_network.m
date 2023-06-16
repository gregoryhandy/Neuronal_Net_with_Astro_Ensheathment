%%
% Creates the connection matrix
% Function Inputs
% Ne: number of excitatory, Ni: number of inhibitory
% Kab: # of connections between neurons of type a and b
%
% Function Outputs
% Vectors of the target/firing neurons pairs and the number of connection
% Connections are chosen with replacement, so neurons can connect
% multiple times
%%
function [target_neuron,firing_neuron,id_start_firing,...
    synapse_num_default,synapse_type] =...
    create_nonSpatial_network(Ne, Ni, Kee, Kie, Kei, Kii)

N = Ne+Ni;

total_num_synapses = Ne*(Kee+Kie)+Ni*(Kei+Kii);
firing_neuron = zeros(total_num_synapses+1,1);
target_neuron = zeros(total_num_synapses,1);
synapse_type = zeros(total_num_synapses,1);
id_start_firing = zeros(Ne+Ni,1);
synapse_num_default = [1:total_num_synapses]';

total_count = 1;
% loop over all neurons
for i = 1:N
    inner_count = 1;
    
    % if the neuron is excitatory 
    if i <= Ne
        
        % excitatory neurons have (Kee+Kie) outgoing connections
        while(inner_count <= (Kee+Kie))
            firing_neuron(total_count) = i;
            
            % choose a random excitatory connection #1-Ne
            if inner_count <= Kee
                target_neuron(total_count) = floor(rand()*Ne)+1;
                synapse_type(total_count) = 1;
            % choose a random inhibitory connection #Ne+1 - Ne+Ni   
            else
                target_neuron(total_count) = Ne+floor(rand()*Ni)+1;
                synapse_type(total_count) = 2;
            end
            
            inner_count = inner_count + 1;
            total_count=total_count+1;
        end
    % if the neuron is inhibitory     
    else
        
        % inhibitory neurons have (Kei+Kii) outgoing connections
        while(inner_count <= (Kei+Kii))
            firing_neuron(total_count) = i;
            
            % choose a random excitatory connection #1-Ne
            if inner_count <= Kei
                target_neuron(total_count) = floor(rand()*Ne)+1;
                synapse_type(total_count) = 3;
            % choose a random inhibitory connection #Ne+1 - Ne+Ni     
            else
                target_neuron(total_count) = Ne+floor(rand()*Ni)+1;
                synapse_type(total_count) = 4;
            end
            
            inner_count = inner_count + 1;
            total_count=total_count+1;
        end
    end
    
    % keep track of the starting index
    % same as, e.g., find(firing_neuron==4320,1)
    if i <= (Ne+1)
        id_start_firing(i) = 1+(Kee+Kie)*(i-1);
    else
        id_start_firing(i) = (Kee+Kie)*(Ne)+1+(Kei+Kii)*(i-Ne-1);
    end

end
    
end

