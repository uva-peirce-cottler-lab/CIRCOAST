% Examine circle packing for high pixel resolution cells and plot results.


output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_ImgDim'];
if isempty(dir(output_path)); mkdir(output_path); end

% Define parameters
umppix = 424.5/512;
img_dim = [64,128, 256, 512,768,1024];
num_trials=1000;
cell_diam_pix = 9;

% Calculate packing ratio and overlap
pack_rat = zeros(numel(cell_diam_pix), num_trials);
for n=1:numel(img_dim)
    n
    trials_tbl = MCM_CirclePacking([img_dim(n) img_dim(n)], num_trials, cell_diam_pix,output_path);
    pack_rat(n,:)=trials_tbl.pack_rat';
end
% trials_tbl.cell_diam_pix = cell_diams_pix';

% Compile Results into table
results_tbl=table();
results_tbl.img_dim  = img_dim';
results_tbl.mean_pack_rat = mean(pack_rat,2);
results_tbl.std_pack_rat = std(pack_rat,[],2);
writetable(results_tbl,[output_path '/hyg_pack_rat.csv']);
% Export results
save([output_path '/MCMRPwoR_Results.mat']);



% Load data, run from here if data is saved to disk
clear all
output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_ImgDim'];
load([output_path '/MCMRPwoR_Results.mat']);



% Generate labels for data
g=bsxfun(@times, img_dim', ones(numel(img_dim),size(pack_rat,2)));

% Test for equal variance across groups
[p,stats] = vartestn(pack_rat(:), g(:));
% Test for
[p_kw,~,stats] = kruskalwallis(pack_rat(:), g(:),...
    arrayfun(@num2str, img_dim, 'Uniform', false));
xlabel('Image Dim. (pix)')
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
set(gcf,'position',[488  245  800  350]);


[c,~,~,gnames] = multcompare(stats);


% Plot mean across groups
u_pack_rat = mean(pack_rat,2);
fprintf('Mean: %.5e, STD %.5e',mean(u_pack_rat), std(u_pack_rat))

% bar(mean(u_pack_rat)-u_pack_rat,mean(u_pack_rat))
