%%
% This code runs a demo for the DiRT simulations used in Handy and 
% Borisyuk, 2023 (corresponding to Fig 1)
%   This code interfaces with mex for speed
%   Also makes use of parfor across trials (see inner loop)
%   Runs over specified values of phi (i.e., fraction of the cleft blocked 
%       by protruding astrocyte)
%   All other default parameters can be adjusted in params_DiRT_2Dcleft
%%
clear; close all; clc;

% Add that folder plus all subfolders to the path.
restoredefaultpath;
folder = fileparts(which('runDiRTSimulation.m')); 
addpath(genpath(folder));
rmpath(folder)

%% Compile the mex code (only needs to happen once)
fprintf('Compiling mex code \n');
mex ./mexCode/DiRT_NeuroAstro_mex.c

%% Specified the values of phi for the simulations
phi_vec = [0 0.4]; % fraction of the cleft blocked by protruding astrocyte

%% Load default parameters and preallocate memory 
params = params_DiRT_2Dcleft(0);
DiRT_results = zeros(params.max_time_points+1,2,params.num_trials,length(phi_vec));

%% Loop over all phi values
for ii = 1:length(phi_vec)
    
    % Load the parameters for this trial
    params = params_DiRT_2Dcleft(phi_vec(ii));
    
    % Seed for the random number generator in the mex code
    randSeed = 833141*ii;
    
    %% Loop over trials 
    % If unable to run in parallel, simply change this to a for-loop
    str1 = sprintf(' %c = %.3f', 981,phi_vec(ii));
    fprintf(strcat('Running the DiRT trials for',str1));
    fprintf(' (%d out of %d values)\n',ii, length(phi_vec));
    output = zeros(params.max_time_points,4,params.num_trials);
    %%
    tic;
    parfor hh = 1:params.num_trials
        % Run the simulation
        [output(:,:,hh)] = DiRT_NeuroAstro_mex(params.max_time,params.dt, params.D, ...
            params.domain_x_min,params.domain_x_max,params.cleft_x_min,params.cleft_x_max,...
            params.domain_y_min,params.domain_y_max,params.cleft_y_min,params.cleft_y_max,...
            params.capture_region_min,params.capture_region_max,...
            params.x_start_loc,params.y_start_loc, params.N_NT, params.tau_r, ...
            params.N_rec, params.receptors, params.print_num, params.max_time_points,...
            params.affinity,randSeed+hh);
    end
    toc;
    
    % Clean up the output a little (can't be in the parfor loop)
    % Some trials might end before max_time is reached (i.e. all particles have 
    % escaped and all receptors recharged). This repeats the end state for the
    % remaining time bins
    for hh = 1:params.num_trials
        temp_index = find(output(:,1,hh)==0,1);
        if ~isempty(temp_index)
            remaining_slots = length(params.t)-temp_index+1;
            output(temp_index:end,2:end,hh) = ...
                repmat(output(temp_index-1,2:end,hh),remaining_slots,1);
        end
        
        % fills in the time column (completely)
        output(:,1,hh) = params.t;
    end
    
    % Save the core results (time and num. of act. receptors)
    % Also concatenate the time = 0 time point
    DiRT_results(:,1,:,ii) = [zeros(1, 1,params.num_trials); output(:,1,:)];
    DiRT_results(:,2,:,ii) = [zeros(1, 1,params.num_trials); params.N_rec - output(:,4,:)];
    
end

%% Ave the results and get curve statistics

max_val = zeros(length(phi_vec),1);
half_max_start= zeros(length(phi_vec),1);
half_max_end= zeros(length(phi_vec),1);
half_max_width= zeros(length(phi_vec),1);
AUC = zeros(length(phi_vec),1);
legend_names = cell(length(phi_vec),1);

for ii = 1:length(phi_vec)
    
    % Average across trials
    DiRT_results_ave = squeeze(mean(DiRT_results,3));
    
    % Half-max width
    [max_val(ii), index_start ] = max(DiRT_results_ave(:,2,ii));
    half_max_start(ii) = find(DiRT_results_ave(:,2,ii) > max_val(ii)/2,1);
    half_max_end(ii) = find(DiRT_results_ave(:,2,ii) < max_val(ii)/2 & ...
        DiRT_results_ave(:,1,ii)>DiRT_results_ave(half_max_start(ii),1,ii),1);
    half_max_width(ii) = DiRT_results_ave(half_max_end(ii),1,ii)-DiRT_results_ave(half_max_start(ii),1,ii);

    % Area under the curve (i.e., synaptic strength)
    AUC(ii) = trapz(DiRT_results_ave(:,1,ii),DiRT_results_ave(:,2,ii));

    legend_names{ii} = sprintf('%.1f',phi_vec(ii));  
end

%% Plot the results
f = figure(113); clf; hold on; h=[];
for ii = 1:length(phi_vec)
    
    h(ii) = plot(DiRT_results_ave(:,1,ii),...
        DiRT_results_ave(:,2,ii),'linewidth',3);
    color_scheme = get(h(ii), 'Color');
    plot(DiRT_results_ave(half_max_start(ii),1,ii), ...
        max_val(ii)/2,'.','markersize',20,'color',color_scheme)
    plot(DiRT_results_ave(half_max_end(ii),1,ii), ...
        max_val(ii)/2,'.','markersize',20,'color',color_scheme)
     
    t_width = [DiRT_results_ave(half_max_start(ii),1,ii)...
        :0.01:DiRT_results_ave(half_max_end(ii),1,ii)];
    plot(t_width,max_val(ii)/2*ones(length(t_width),1),'--',...
        'linewidth',2,'color',color_scheme)
    
    str1 = sprintf('For %c = %.2f, the strength of the synapse is %.2f', 981,phi_vec(ii),AUC(ii));
    str2 = sprintf(' and the half-width max is %.2f', half_max_width(ii));
    fprintf(strcat(str1,str2)); fprintf('\n');
end
set(gca,'fontsize',16)
xlim([0 1])
legend(h,strcat('\phi= ',legend_names))
xlabel('Time (a.u.)')
ylabel('# of Activated Receptors')
yticks([0:10:50])

%% Save the results
save('./SimResults/DiRT_Results','DiRT_results','params','phi_vec','-v7.3')
