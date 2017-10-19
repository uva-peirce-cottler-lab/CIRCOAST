
clear all
clc
vessel_diam_um = 9;
cell_diam_um = 15;
umppix = 1.26;
img_dim = [1024 1024];
max_trials_per_cycle = 300;

vessel_densities = [.005 .03 .05 .1 .2 .3 .4 .5 .6  .7 .8 .9 .95];
events_per_trial_vector = [1 2 3 5 50 100 500 1000 2500 7500 10000 15000 20000 25000 30000 60000];
tot_trials = 10000;
handles = [];


proj_path = getappdata(0,'proj_path');
output_path = [proj_path '/output'];
mkdir(output_path);

output_name = 'benchmark_output.csv'; 

c = 1; figure('units', 'normalized');
ha = tight_subplot(numel(vessel_densities), numel(events_per_trial_vector), 0, 0, 0);
for n = 1:numel(vessel_densities)
    
    fprintf('Generating Vessel Network: %.f / %.f\n', n, numel(vessel_densities));
    bw_vessels = VesselGenerator_SpawnHorizontal(img_dim , vessel_diam_um, umppix, ...
        'VesselDensityFraction', vessel_densities(n));
    
    for k = 1:numel(events_per_trial_vector)
        
        % Estimate # trials that can be computed concurrently in memory
        [comp_img, max_trials_per_cycle, ~] = ArcasGui_RunSingleTrial(bw_vessels,...
            events_per_trial_vector(end),1,cell_diam_um,umppix);
        
%       % Run simulation
        fprintf('Running Simulation: %.f / %.f\n', k, numel(events_per_trial_vector));
        [sim_frac_means(n,k), sim_frac_stds(n,k), event_vasc_dists, comp_img,...
            ~, mean_frac_cell_density(n,k), mean_frac_cell_overlaps(n,k)] = ...
            ArcasGui_monteCarloSim_Driver(bw_vessels, cell_diam_um, umppix, ...
            tot_trials, events_per_trial_vector(k), max_trials_per_cycle, handles);
        
        
       out_st = ...
           CELLCOAV_BMRP(bw_vessels, cell_diam_um, umppix,events_per_trial_vector(k));
        binom_frac_means(n,k) = out_st.binom_frac_mean;
        binom_frac_stds(n,k) = out_st.binom_frac_std;
%        
%        out_st = ...
%     CELLCOAV_HMRP(bw_vessels, cell_diam_um, umppix, n_obs_success,events_per_trial_vector(k))


       
        img_name = ['img vd_' num2str(vessel_densities(n)) ...
            ' events_' num2str(events_per_trial_vector(k)) ...
            ' eventdiam_' num2str(cell_diam_um) '.tif'];
        
        
        subaxis(numel(vessel_densities),numel(events_per_trial_vector),...
            c, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
        imshow(comp_img);
        c = c + 1;
        
        imwrite(comp_img, [output_path '/' img_name]);
    end
    
end
% fclose(fid);
% Save all variables to disk
save('benchmark_data.mat');
load('benchmark_data.mat');
% keyboard


% Recreate Plot figure
c = 1; figure('units', 'normalized');
ha = tight_subplot(numel(vessel_densities), numel(events_per_trial_vector), 0, 0, 0);
for n = 1:numel(vessel_densities)
    fprintf('Generating Vessel Network: %.f / %.f\n', n, numel(vessel_densities));
    bw_vessels = VesselGenerator_SpawnHorizontal([1024 1024], vessel_diam_um, umppix, ...
        'VesselDensity', vessel_densities(n));
    for k = 1:numel(events_per_trial_vector)
        [~, ~, ~, comp_img] = ...
            ArcasGui_monteCarloSim_Driver(bw_vessels, cell_diam_um, umppix, ...
            1, events_per_trial_vector(k), 1, handles);
        subaxis(numel(vessel_densities),numel(events_per_trial_vector),...
            c, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
        imshow(comp_img); c = c + 1;
    end
end
htile = gcf;
saveas(gcf,'tile.fig')
openfig('tile.fig')


% Compare monte carlo and Binomial distribution, row 6, 1-6
figure;
row =4; ind = 1:9;
errorbar(binom_frac_means(row,ind), binom_frac_stds(row,ind), 'r--')
hold on
errorbar(sim_frac_means(row,ind), sim_frac_stds(row,ind), 'b--')
set(gca,'XTickLabel', horzcat({''}, ...
    strread(num2str(events_per_trial_vector(ind)),'%s')',{''}))
ay=ylim; ax=xlim;
axis([ax 0 ay(2)]); grid off 
ylabel('Mean Colocalization Fraction'); xlabel('Cells per FOV');
legend({'Binomial','Monte Carlo'})
% title(sprintf('Comparison of Model Predictions of Random Colocalization\n from Simulated Images'))
hold off
%% Grab images used for this comparison
imgs={}; figure(htile);
for n = 1:numel(ind)
   subaxis((numel(events_per_trial_vector)*(row-1))+n);
   imgs{n} = getimage(gca);
end
% Create subplot
figure;
for n = 1:numel(ind)
   subaxis(1,numel(ind),n,'Spacing', 0.005, 'Padding', 0, 'Margin', 0);
   imshow(imgs{n},'Parent',gca); 
end



% Compare monte carlo and Binomial distribution, row 6, 1-6
figure;
col =6; ind = 1:9;
errorbar(binom_frac_means(ind,col), binom_frac_stds(ind,col), 'r--')
hold on
errorbar(sim_frac_means(ind,col), sim_frac_stds(ind,col), 'b--')
set(gca,'XTickLabel', horzcat({''}, ...
    strread(num2str(events_per_trial_vector(ind)),'%s')',{''}))
ay=ylim; ax=xlim;
axis([ax 0 ay(2)]); grid off 
ylabel('Mean Colocalization Fraction'); xlabel('Vessel Density (um/mm^2)');
legend({'Binomial','Monte Carlo'})
% title(sprintf('Comparison of Model Predictions of Random Colocalization\n from Simulated Images'))
hold off
%% Grab images used for this comparison
imgs={}; figure(htile);
for n = 1:numel(ind)
   subaxis(col+ind(n)*numel(events_per_trial_vector));
   imgs{n} = getimage(gca);
end
% Create subplot
figure;
for n = 1:numel(ind)
   subaxis(1,numel(ind),n,'Spacing', 0.005, 'Padding', 0, 'Margin', 0);
   imshow(imgs{n},'Parent',gca); 
end






figure; plot(vessel_densities,binom_frac_means(:,1),'k','LineWidth', 3)
hold on; plot(vessel_densities,sim_frac_means); hold off
legend(horzcat({'Binom. Est.'},...
    strcat(regexp(num2str(events_per_trial_vector),'\s*','split'),{' Cells'})))
h = sort(get(gca,'Children'),'ascend');
set(h(9:15),'LineStyle','--'); set(h(16:17),'LineStyle',':');
title('Mean Colocalization Fraction Over Range of Vessel Densities')
ylabel('Mean Colocalization Fraction'); xlabel('Vessel Density Fraction');
% set(gca,'XTickLabel', vessel_densities);
grid on
saveas(gcf, 'coloc_mean.fig');

figure
plot(vessel_densities,sim_frac_stds);
legend(strcat(regexp(num2str(events_per_trial_vector),'\s*','split'),{' Cells'}))
h = sort(get(gca,'Children'),'ascend');
set(h(8:14),'LineStyle','--'); set(h(15:16),'LineStyle',':');
title('Monte Carlo STD of Colocalization Fraction Over Range of Vessel Densities')
ylabel('STD Colocalization Fraction'); xlabel('Vessel Density Fraction');
% set(gca,'XTickLabel', vessel_densities);
grid on
saveas(gcf, 'monte_carlo_coloc_std.fig');


figure
plot(vessel_densities,binom_frac_stds);
legend(strcat(regexp(num2str(events_per_trial_vector),'\s*','split'),{' Cells'}))
h = sort(get(gca,'Children'),'ascend');
set(h(8:14),'LineStyle','--'); set(h(15:16),'LineStyle',':');
title('Binomial STD of Colocalization Fraction Over Range of Vessel Densities')
ylabel('STD Colocalization Fraction'); xlabel('Vessel Density Fraction');
% set(gca,'XTickLabel', vessel_densities); 
grid on
saveas(gcf, 'binom_coloc_std.fig');

figure;
plot(events_per_trial_vector(1:11),mean_frac_cell_density(1,1:11))
title('Mean Cell Density Fraction Over Range of Injected Cells')
ylabel('Mean Cell Density Fraction'); xlabel('Cells per FOV');
% set(gca,'XTickLabel', vessel_densities); 
grid on


figure;
plot(events_per_trial_vector(1:11),mean_frac_cell_overlaps(1,1:11) - 1./events_per_trial_vector(1:11))
title('Fraction of Cells Overlapping Over Range of Cell Count')
ylabel('Mean Fraction Cells Overlaping'); xlabel('Cells per FOV');
% set(gca,'XTickLabel', vessel_densities); 
grid on
% t = 0:0.001:2*pi+0.001;
% figure(2);
% for i = 1 : 25;
%     subaxis(8,7,i, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0);
%     plot(t, sin(i*t));
%     axis tight
%     axis off
% end





