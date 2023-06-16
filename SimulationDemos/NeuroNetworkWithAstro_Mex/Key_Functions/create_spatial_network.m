%%
%  Creates a spatial network of neurons, using the design by
%  Rosenbaum et al., 2017
%%

function [Wex1, Wex2, Wix1, Wix2,target_neuron, firing_neuron,...
          id_start_firing, synapse_num_default,synapse_type]=...
          create_spatial_network(Ne1, Kee, betaee, Kie, betaie,...
                                 Ni1, Kii, betaii, Kei, betaei,...
                                 Nx1, Kex, betaex, Kix, betaix)

% Parameters
Ne = Ne1^2; Ni = Ni1^2; Nx = Nx1^2;


total_num_synapses = Ne*(Kee+Kie)+Ni*(Kei+Kii);
firing_neuron = zeros(total_num_synapses+1,1);
target_neuron = zeros(total_num_synapses,1);
synapse_type = zeros(total_num_synapses,1);
id_start_firing = zeros(Ne+Ni,1);
synapse_num_default = [1:total_num_synapses]';

% Preallocate memory
Wee1 = zeros(Ne*Kee,1); Wee2 = zeros(Ne*Kee,1); Wee = zeros(Ne*Kee,1);
Wie1 = zeros(Ne*Kie,1); Wie2 = zeros(Ne*Kie,1); Wie = zeros(Ne*Kie,1);
Wii1 = zeros(Ni*Kii,1); Wii2 = zeros(Ni*Kii,1); Wii = zeros(Ni*Kii,1);
Wei1 = zeros(Ni*Kei,1); Wei2 = zeros(Ni*Kei,1); Wei = zeros(Ni*Kei,1);

Wex1 = zeros(Nx*Kex,1); Wex2 = zeros(Nx*Kex,1);
Wix1 = zeros(Nx*Kix,1); Wix2 = zeros(Nx*Kix,1);

% Loop through the excitatory neurons
for j=0:(Ne-1)
    % Find index of cell along each dimension */
    j1 = floor(j/Ne1);
    j2 = mod(j,Ne1);
    
    % Generate vectors of exc and inh postsynaptic targets */
    % Inputs to CircRandNfun : mu, sigma, rand_min, rand_max, n)
    % x coordinate of neuron
    Wee1((1+Kee*j):Kee*(j+1)) = CircRandNfun(j1,betaee,0,Ne1-1,Kee);
    % y coordinate of neuron
    Wee2((1+Kee*j):Kee*(j+1)) = CircRandNfun(j2,betaee,0,Ne1-1,Kee);
    %Wei1[(j-Ne)*Kei+k]*Ne1+Wei2[(j-Ne)*Kei+k]
    Wee((1+Kee*j):Kee*(j+1)) = Wee1((1+Kee*j):Kee*(j+1))*Ne1+Wee2((1+Kee*j):Kee*(j+1));
    
    % x coordinate of neuron
    Wie1((1+Kie*j):Kie*(j+1)) = CircRandNfun(j1*(Ni1/Ne1),betaie,0,Ni1-1,Kie);
    %y coordinate of neuron
    Wie2((1+Kie*j):Kie*(j+1))  = CircRandNfun(j2*(Ni1/Ne1),betaie,0,Ni1-1,Kie);
    
    Wie((1+Kie*j):Kie*(j+1)) = Ne+Wie1((1+Kie*j):Kie*(j+1))*Ni1+Wie2((1+Kie*j):Kie*(j+1));
    
    firing_neuron((1+(Kee+Kie)*j):(Kee+Kie)*(j+1)) = j*ones(Kee+Kie,1);
    target_neuron((1+(Kee+Kie)*j):(Kee+Kie)*(j+1)) = [Wee((1+Kee*j):Kee*(j+1));...
        Wie((1+Kie*j):Kie*(j+1))];

    synapse_type((1+(Kee+Kie)*j):(Kee+(Kee+Kie)*j)) = 1;
    synapse_type((1+Kee+(Kee+Kie)*j):(Kee+Kie)*(j+1)) = 2;
end


% loop through the inhibitory neurons
for j=0:(Ni-1)
    
    % Find index of cell along each dimension
    j1 = floor(j/Ni1);
    j2 = mod(j,Ni1);
    
    % Generate vectors of exc and inh postsynaptic targets
    Wei1((1+Kei*j):Kei*(j+1)) = CircRandNfun(j1*(Ne1/Ni1),betaei,0,Ne1-1,Kei);
    Wei2((1+Kei*j):Kei*(j+1)) = CircRandNfun(j2*(Ne1/Ni1),betaei,0,Ne1-1,Kei);
    
    Wei((1+Kei*j):Kei*(j+1)) =  Wei1((1+Kei*j):Kei*(j+1))*Ne1+Wei2((1+Kei*j):Kei*(j+1));
    
    Wii1((1+Kii*j):Kii*(j+1)) = CircRandNfun(j1,betaii,0,Ni1-1,Kii);
    Wii2((1+Kii*j):Kii*(j+1)) = CircRandNfun(j2,betaii,0,Ni1-1,Kii);
    
    Wii((1+Kii*j):Kii*(j+1)) = Ne+Wii1((1+Kii*j):Kii*(j+1))*Ni1+Wii2((1+Kii*j):Kii*(j+1));
       
    firing_neuron(((Kee+Kie)*Ne+1+(Kei+Kii)*j):((Kee+Kie)*Ne+(Kei+Kii)*(j+1)))...
        = Ne+j*ones(Kei+Kii,1);
    target_neuron(((Kee+Kie)*Ne+1+(Kei+Kii)*j):((Kee+Kie)*Ne+(Kei+Kii)*(j+1)))...
        = [Wei((1+Kei*j):Kei*(j+1)); Wii((1+Kii*j):Kii*(j+1))];
    
    synapse_type(((Kee+Kie)*Ne+1+(Kei+Kii)*j):((Kee+Kie)*Ne+Kei+(Kei+Kii)*j)) = 3;
    synapse_type(((Kee+Kie)*Ne+1+Kee+(Kei+Kii)*j):((Kee+Kie)*Ne+(Kei+Kii)*(j+1))) = 4;
end


% keep track of the starting index
% same as, e.g., find(firing_neuron==4320,1)
for i = 0:(Ne+Ni-1)    
    if i < Ne
        id_start_firing(i+1) = (Kee+Kie)*i;
    else
        id_start_firing(i+1) = (Kee+Kie)*(Ne)+(Kei+Kii)*(i-Ne);
    end
end

% Keep the indexing the same as Rosenbaum et al. for the feedforward
% connections (might need to change)
% loop through the feedforward connections
for j=0:(Nx-1)
    
    % Find index of cell along each dimension
    j1= floor(j/Nx1);
    j2= mod(j,Nx1);
    
    % Generate vectors of exc and inh postsynaptic targets
    Wex1((1+Kex*j):Kex*(j+1)) = CircRandNfun(j1*(Ne1/Nx1),betaex,0,Ne1-1,Kex);
    Wex2((1+Kex*j):Kex*(j+1)) = CircRandNfun(j2*(Ne1/Nx1),betaex,0,Ne1-1,Kex);
    
    Wix1((1+Kix*j):Kix*(j+1)) = CircRandNfun(j1*(Ni1/Nx1),betaix,0,Ni1-1,Kix);
    Wix2((1+Kix*j):Kix*(j+1)) = CircRandNfun(j2*(Ni1/Nx1),betaix,0,Ni1-1,Kix);
end

end