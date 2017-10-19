function [comp_img, max_trials_per_cycle, minutes_for_simulation, tot_exec_bytes] = ...
    ArcasGui_RunSingleTrial(bw_img,events_per_trial,tot_trials,event_diam,umppix)
% pixel diameter of agent
event_pix_diam = ceil(event_diam/umppix);

% Estimate time for all trials by running one trial (overestimated)
tic;
[~, ~, ~,comp_img,frac_hits] = CELLCOAV_MCMRP(bw_img, 1, ...
    events_per_trial, event_pix_diam);
timerVal = toc;
% frac_hits*events_per_trial
% keyboard
% imwrite(comp_img,'a.tif')
% fprintf('frac_hits: %.3f\n',frac_hits)

%% Estimate Memory
% Matlab profiler fails to proprly estimate memory for a given mfile,
% leaves out many function calls such as randi() for no reason.
trial_vect = [1 100];
for n = 1:2
    % % Estimate memory for one trial
    c = onCleanup(@(x) profile('off'));
    % Run profile viewer while executing one trial, get mem usage
    afun = @(x,y,z,e) CELLCOAV_MCMRP(x,y,z,e);
    profile -history -memory on
    [~, ~, ~,~, ~, var_mem(n)] = afun(bw_img, trial_vect(n), events_per_trial, event_pix_diam);
    profile off
    prof_text = profile('info');
    % Calculate memory required to run one trial
    tot_exec_time(n) = max(arrayfun(@(x) x.TotalTime, prof_text.FunctionTable));
    tot_exec_bytes(n) = max(arrayfun(@(x) x.PeakMem, prof_text.FunctionTable));
end
 
% Query available memory in matlab
[userview] = memory;
total_mem = tot_exec_bytes + var_mem;
% Have to xshift memory vals beacuse of badly condiitoned input for polyfit
p = polyfit(total_mem-min(total_mem),trial_vect,1);
tot_trials = polyval(p,userview.MemAvailableAllArrays+min(total_mem));
% Calculate max number of trials that can be run in matlab without exceeding
% System memory (with safety% maargin)
max_trials_per_cycle = tot_trials * 0.1;

% Calculate Time with fitting line to 2 test runs
p = polyfit(trial_vect,tot_exec_time,1);
minutes_for_simulation = abs(polyval(p,tot_trials)/60);






