function [sim_hit_mean, sim_hit_std, event_vasc_dists, comp_img, trials_hit_frac,...
    mean_frac_cell_density, mean_frac_cell_overlaps] = ...
    ArcasGui_monteCarloSim_Driver(img, event_diam_um, umppix, ...
    tot_trials, events_per_trial, max_trials_per_cycle, handles)

if ~exist('handles','var'); handles = []; end

% pixel diameter of agent
event_pix_diam = event_diam_um/umppix;

% Run a series of simulations with only a portion of trails completed
trials_per_cycle = ceil(max_trials_per_cycle/1.5);

% Initialize status bar
if ~isempty(handles) 
%     sb = statusbar('  Starting Simulation');
%     set(sb.ProgressBar, 'Minimum',0, 'Maximum',100, 'Value',0);
%     sb.ProgressBar.setVisible(true);
end
% To run the specified number of simulations, the trials are split into
% groups that can fit into the local computer's RAM
for n = 1:ceil(tot_trials/trials_per_cycle)
    
    % Check if cancel button was pressed
    if ~isempty(handles) && get(handles.cancel_togglebutton,'Value') == 1;
        [sim_hit_mean, sim_hit_std, event_vasc_dists, comp_img] = deal([]);
        statusbar;
        return;
    end
    
%     fprintf('Cycle %.f / %.f \n', n, ceil(tot_trials/trials_per_cycle));
    if n == ceil(tot_trials/trials_per_cycle);
        ntrials = tot_trials - (n-1)*trials_per_cycle;
    else
        ntrials =  trials_per_cycle;
    end
    
    % Run sub simulation
%     [~, ~, event_vasc_dists_cell{n}, comp_img, trials_hit_frac_cell{n}, ...
%         frac_cell_densities_cell{n}, frac_cell_overlaps_cell{n}] = ...
%         ArcasGui_monteCarloSim_MedMem(img, ntrials, events_per_trial, event_pix_diam);
    
    % Run sub simulation
    [~, ~, event_vasc_dists_cell{n}, comp_img, trials_hit_frac_cell{n}] = ...
        CELLCOAV_MCMRP(img, ntrials, events_per_trial, event_pix_diam);
    
    
    if ~isempty(handles)
        imshow(comp_img,'Parent', handles.img_axes);
%         statusbar(handles.ArcasGui_figure, '  %.f%%:  %.f Trials Complete of %.f', ...
%             n*trials_per_cycle/tot_trials*100, n*trials_per_cycle, tot_trials)
%         sb.ProgressBar.setValue(n/ceil(tot_trials/trials_per_cycle) * 100); drawnow; pause(.01);
    end
end
% Delete Status Bar
if ~isempty(handles)
%     statusbar;
end

trials_hit_frac = vertcat(trials_hit_frac_cell{:});
% mean_frac_cell_density  = mean(vertcat(frac_cell_densities_cell{:}));
% mean_frac_cell_overlaps = mean(vertcat(frac_cell_overlaps_cell{:}));
mean_frac_cell_density  = 1;
mean_frac_cell_overlaps = 1;

% Simulation mean and std for hits per trial
% TODO calculate std from subgroups with known subgroup mean and std
% https://www.youtube.com/watch?v=xP4I81FL9jA
sim_hit_mean = mean(trials_hit_frac);
sim_hit_std = std(trials_hit_frac);

% keyboard
event_vasc_dists = horzcat(event_vasc_dists_cell{:});


