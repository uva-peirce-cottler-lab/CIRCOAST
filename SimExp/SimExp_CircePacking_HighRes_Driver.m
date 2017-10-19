% Examine circle packing for high pixel resolution cells and plot results.


output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_HighRes'];
if isempty(dir(output_path)); mkdir(output_path); end

% Define parameters
umppix = 424.5/512;
img_dim = [1024 1024];
num_trials=1000;
cell_diams_pix = [55:6:101 141:10:191 201 225 251 301 351 401];
cell_diams_um = round(cell_diams_pix.* umppix*10)/10;

% Calculate packing ratio and overlap
pack_rat = zeros(numel(cell_diams_pix), num_trials);
for n=1:numel(cell_diams_pix)
    n
    trials_tbl = MCM_CircePacking(img_dim, num_trials, cell_diams_pix(n),output_path);
    pack_rat(n,:)=trials_tbl.pack_rat';
end
% trials_tbl.cell_diam_pix = cell_diams_pix';

% Compile Results into table
hyg_pac_rat_tbl=table();
hyg_pac_rat_tbl.cell_diam_pix = cell_diams_pix';
hyg_pac_rat_tbl.mean_pack_rat = mean(pack_rat,2);
hyg_pac_rat_tbl.std_pack_rat = std(pack_rat,[],2);
writetable(hyg_pac_rat_tbl,[output_path '/hyg_pack_rat.csv']);
save([output_path '/MCMRPwoR_Results.mat']);



% Load data, run from here if data is saved to disk
clear all
output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest_HighRes'];
load([output_path '/MCMRPwoR_Results.mat']);


% If required, only examine section of data
% pack_rat=pack_rat(1:16,:);

% Generate labels for data
lbl=bsxfun(@times, cell_diams_um', ones(numel(cell_diams_pix),size(pack_rat,2)));

% Test for equal variance across groups
[p,stats] = vartestn(pack_rat(:), lbl(:));
% Test for
[p_kw,~,stats] = kruskalwallis(pack_rat(:), lbl(:),...
    arrayfun(@num2str, cell_diams_um, 'Uniform', false));
xlabel('Cell Diam. (um)')
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




%% CDVF*p approximationg of k Test
vect_vld= [30 25 20 17.5 15 12.5 10 7.5 5 2.5];

results_tbl=table();
results_tbl.img_ind=(1:50)';
results_tbl.cdvf= zeros(50,1);
results_tbl.coloc_cells = zeros(50,1);
results_tbl.max_cells = zeros(50,1);
results_tbl.vld = zeros(50,1);
results_tbl.vaf = zeros(50,1);




k=1;
for r=6:10
    % Generate images
    [~, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
        3*6, 'VesselLengthDensityLimits', [40 0]);
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel);
    
    % Iteratively remove vessel segments to a range of vessel densities
    [ves_cell{k}, skel_cell{k}, stat_st] =...
        VesselGen_RegressNetwork(bw_skel, 3,umppix, ...
        'VesselLengthDensityTarget',vect_vld(1),'BW_Protected',bw_protect);
    results_tbl.vld(k) = stat_st.vld_mmpmm2;
    results_tbl.vaf(k) = stat_st.vf;
    for c=2:numel(vect_vld)
        [ves_cell{k+c-1}, skel_cell{k+c-1}, stat_st] =...
            VesselGen_RegressNetwork(skel_cell{k+c-2}, 3,umppix, ...
            'VesselLengthDensityTarget',vect_vld(c),'BW_Protected',bw_protect);
        results_tbl.vld(k+c-1) = stat_st.vld_mmpmm2;
    results_tbl.vaf(k+c-1) = stat_st.vf;
    end
    k=k+numel(vect_vld);
end

% Write images to disk
for n=1:numel(ves_cell)
    imwrite(ves_cell{n}, [output_path '/bw_vessels.tif'],'WriteMode','append');
end

% Load vascular Images
img_info=imfinfo([output_path '/bw_vessels.tif']);
bw_vessels=false(img_info(1).Height,img_info(1).Width,length(img_info));
for z=1:numel(img_info)
    bw_vessels(:,:,z)=imread([output_path '/bw_vessels.tif'],'Index',z,...
        'Info',img_info);
end


% 
% % Calculate CDVF for each image
% for z=1:size(bw_vessels,3)
%    cdvf(z) =  ArcasGui_PercentCellOverlap(bw_vessels(:,:,z),cell_diam_um,umppix)
% end


% For each cell size
%          For each vascular network image: average difference between
%                  cdvf*N          coloc_cells
%     across 100 images 
% cell_size, vasc image
avg_cell_ov_diff = zeros(numel(cell_diams_pix),size(bw_vessels,3));

for c=1:numel(cell_diams_pix)
    fprintf('%i',c);
    for z=1:size(bw_vessels,3)
        fprintf('.');
        % Load all cell images for this cell size
        bw_cctr_trials = Load_CellPlacement_Series(['bw_cell_CD' num2str(cell_diams_pix(c)) '*'],...
            output_path,1);
        
       % CDVF 1xnum_trials vector
       cdvf = ArcasGui_PercentCellOverlap(bw_vessels(:,:,z),cell_diams_pix(c),1);
       % Get total cells found in each image (sum # cell centers)
       N = squeeze(sum(sum(bw_cctr_trials,1),2));
       num_coloc_cells_cdvf = floor(N.*cdvf);
       
       % For each image, find distance from ctr point to near vessel, calc
       % fraction of colocalizaing cells
       ed_vessel =bwdist(bw_vessels(:,:,z));
       ed_ctr_dist=bsxfun(@times, ed_vessel,single(bw_cctr_trials));
       coloc_ctr = (ed_ctr_dist <= (cell_diams_pix(c)-1)/2) & (ed_ctr_dist>0);
       
       num_coloc_cells_mcm = squeeze(sum(sum(coloc_ctr,1),2));
       
       avg_cell_ov_diff (c,z) = mean(num_coloc_cells_cdvf- num_coloc_cells_mcm);
    end
    fprintf('\n');
end

   % Generate labels for data
lbl=bsxfun(@times, cell_diams_um', ones(numel(cell_diams_pix),size(bw_vessels,3)));

