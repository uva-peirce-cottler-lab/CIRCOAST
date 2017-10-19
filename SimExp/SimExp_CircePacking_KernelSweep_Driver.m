% Examine circle packing for high pixel resolution cells and plot results.


output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_KernelSweep'];
if isempty(dir(output_path)); mkdir(output_path); end

% Define parameters
umppix = 424.5/512;
img_dim = [512 512];
num_trials=1000;
cell_diams_pix = [3:2:55];
cell_diams_um = round(cell_diams_pix.* umppix*10)/10;

% Calculate packing ratio and overlap
pack_rat = zeros(numel(cell_diams_pix), num_trials);
for n=1:numel(cell_diams_pix)
    n
    trials_tbl = MCM_CirclePacking(img_dim, num_trials, cell_diams_pix(n),output_path);
    pack_rat(n,:)=trials_tbl.pack_rat';
end
% trials_tbl.cell_diam_pix = cell_diams_pix';

% Compile Results into table
results_tbl=table();
results_tbl.cell_diam_pix = cell_diams_pix';
results_tbl.mean_pack_rat = mean(pack_rat,2);
results_tbl.std_pack_rat = std(pack_rat,[],2);
writetable(results_tbl,[output_path '/hyg_pack_rat.csv']);
% Export results
save([output_path '/MCMRPwoR_Results.mat']);



% Load data, run from here if data is saved to disk
clear all
output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_KernelSweep'];
load([output_path '/MCMRPwoR_Results.mat']);


% If required, only examine section of data
pack_rat=pack_rat(1:17,:);
cell_diams_pix = cell_diams_pix(1:17);

% Generate labels for data
g=bsxfun(@times, cell_diams_pix', ones(numel(cell_diams_pix),size(pack_rat,2)));

% Test for equal variance across groups
[p,stats] = vartestn(pack_rat(:), g(:));
% Test for
[p_kw,~,stats] = kruskalwallis(pack_rat(:), g(:),...
    arrayfun(@num2str, cell_diams_pix, 'Uniform', false));
xlabel('Cell Diam. (pix)')
ylabel('Packing Ratio')

hTitle=get(gca,'title');
hLegend = get(gca,'Legend');
hXLabel = get(gca,'XLabel');
hYLabel = get(gca, 'YLabel');
set(gca                       , ...
    'FontName'   , 'Arial' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Ariel');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 14          , ...
    'FontWeight' , 'bold'      );
set(gcf,'position',[488  245  850  350]);


[c,~,~,gnames] = multcompare(stats);


% Plot mean across groups
u_pack_rat = mean(pack_rat,2);
fprintf('Mean: %.5e, STD %.5e',mean(u_pack_rat), std(u_pack_rat))

% bar(mean(u_pack_rat)-u_pack_rat,mean(u_pack_rat))
