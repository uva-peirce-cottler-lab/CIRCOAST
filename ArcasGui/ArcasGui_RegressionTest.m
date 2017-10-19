
% Test to make sure you get similiar values with changes to monte carlo
% program


% clear all
% clc
vessel_diam_um = 9;
cell_diam_um = 15;
umppix = 1.26;
img_dim = [1024 1024];
max_trials_per_cycle = 300;

vessel_densities = [.005 .03 .05 .1 .2 .3 .4 .5 .6  .7 .8 .9 .95];
events_per_trial_vector = [1 2 3 5 50 100 500 1000 2500 7500 10000 15000 20000 250000 30000 60000];
tot_trials = 1000;
handles = [];



event_pix_diam = ceil(cell_diam_um/umppix);



% fprintf('Generating Vessel Network: %.f / %.f\n', n, numel(vessel_densities));
bw_vessels = VesselGenerator_SpawnHorizontal(img_dim , vessel_diam_um, umppix, ...
    'VesselDensityFraction', .2);



[hit_mean, hit_std,~,~,trials_hit_frac_cell] =...
    ArcasGui_monteCarloSim_MedMem(bw_vessels, 1000, 1, event_pix_diam)



