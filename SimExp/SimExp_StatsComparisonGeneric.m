function [] = SimExp_StatsComparisonGeneric()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ClearPreviousData=1;
%Output Directory
output_path = [getappdata(0,'proj_path') '/temp_data/StatsComparison'];
mkdir(output_path);
delete([output_path '/*.*']);

% Simulation Parameters
N_tot=10;
n_tot=4;
grp_tot = 2;
img_dim = [512 512];
umppix = 424.5/512;
cell_diam_um = 15;
cell_diam_pix = cell_diam_um/umppix;


% Data table
if isempty(dir([output_path '/arcas_results.csv']))
    % image_Name, group_id, sample_id, rep_id, vld, tot_cells, obs_col_cells, obs_icf, bmrp_icf
    tbl=table;
    tbl.img_name = cell(grp_tot*N_tot*n_tot,1);
    tbl.group_id = vertcat(ones(N_tot*n_tot,1), 2.*ones(N_tot*n_tot,1));
    tbl.sample_id = zeros(grp_tot*N_tot*n_tot,1);
    tbl.repl_id = zeros(grp_tot*N_tot*n_tot,1);
    tbl.vld = zeros(grp_tot*N_tot*n_tot,1);
    tbl.tot_cells = round(normrnd(vertcat(ones(N_tot*n_tot,1).*30, ones(N_tot*n_tot,1).*25),2));
    tbl.obs_icf = zeros(grp_tot*N_tot*n_tot,1);
    tbl.bmrp_icf = zeros(grp_tot*N_tot*n_tot,1);
    tbl.bmrp_cellcoav_p = zeros(grp_tot*N_tot*n_tot,1);
    tbl.hmrp_icf = zeros(grp_tot*N_tot*n_tot,1);
    tbl.hmrp_cellcoav_p = zeros(grp_tot*N_tot*n_tot,1);
    
    % Requested vld is not exactly the same as actual
    req_vld = normrnd(vertcat(ones(N_tot*n_tot,1).*30, ones(N_tot*n_tot,1).*20),2);
else
    tbl=readtable([output_path '/arcas_results.csv']);
end

% Create image dataset
hw = waitbar(0,'Processing Images');
k=1;
for N=1:N_tot
    for n=1:n_tot
        fprintf('N%d_n%d.tif',N,n);
%         keyboard
        % If high density image does not exist, create seed vessel network
        if isempty(tbl.img_name{k})
            % Create Initial over dense vessel network
            [~, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
                3*6, 'VesselLengthDensityLimits', [50 0],'FigureVisible','off');
            % Identify central vessel to protect
            bw_protect = VesselGen_SelectProtectedSegment(bw_skel,'FigureVisible','off');
        end
        
        % *** Group 1
        ind = k;
        % Iteratively remove vessel segments to a range of vessel densities
        if isempty(tbl.img_name{ind})
            [bw_vessel, ~, stat_st] =...
                VesselGen_RegressNetwork(bw_skel, 3,umppix, ...
                'VesselLengthDensityTarget',req_vld(ind),...
                'BW_Protected',bw_protect,'FigureVisible','off');
            
            tbl.vld(ind) = stat_st.vld_mmpmm2;
            % Populate cells in image
            [~, ~, ~,comp_img,tbl.obs_icf(ind)] = CELLCOAV_MCMRP(bw_vessel, 1, ...
                tbl.tot_cells(ind), cell_diam_pix);
            % Export Image
            tbl.img_name{ind} =  sprintf('G1_N%d_n%d.tif',N,n);
            imwrite(comp_img,[output_path '/' tbl.img_name{ind}]);
        else
            comp_img = imread([output_path '/' tbl.img_name{ind}]);
            bw_vessel=comp_img(:,:,2);
        end
        % Get binomial stats
        binom_st = CELLCOAV_BMRP(bw_vessel, ...
            cell_diam_um, umppix,tbl.tot_cells(ind),'n_obs_success',...
            tbl.obs_icf(ind)*tbl.tot_cells(ind));
        tbl.bmrp_icf(ind) = binom_st.binom_frac_mean;
        tbl.bmrp_cellcoav_p(ind) = binom_st.binom_p_nerandom;
        % Hypergeometric model of random cell placement
        hyg_st = CELLCOAV_HMRP(bw_vessel, ...
            cell_diam_um, umppix,tbl.tot_cells(ind),'n_obs_success',...
            tbl.obs_icf(ind)*tbl.tot_cells(ind));
        tbl.hmrp_icf_mean(ind) = hyg_st.hmrp_icf_mean;
        tbl.hmrp_cellcoav_p(ind) = hyg_st.hmrp_cellcoav_p;
        % Fill in other metadata
        tbl.sample_id(ind) = N;
        tbl.repl_id(ind) = n;
        
        
        
        % *** Group 2
        ind = k+N_tot*n_tot;
        % Iteratively remove vessel segments to a range of vessel densities
        if isempty(tbl.img_name{ind})
            [bw_vessel, ~, stat_st] =...
                VesselGen_RegressNetwork(bw_skel, 3,umppix, ...
                'VesselLengthDensityTarget',req_vld(ind),...
                'BW_Protected',bw_protect,'FigureVisible','off');
            tbl.vld(ind) = stat_st.vld_mmpmm2;
            % Populate cells in image
            [~, ~, ~,comp_img,tbl.obs_icf(ind)] = CELLCOAV_MCMRP(bw_vessel, 1, ...
                tbl.tot_cells(ind), cell_diam_pix);
            % Export Image
            tbl.img_name{ind} =  sprintf('G2_N%d_n%d.tif',N,n);
            imwrite(comp_img,[output_path '/' tbl.img_name{ind}]);
        else
            comp_img = imread([output_path '/' tbl.img_name{ind}]);
            bw_vessel=comp_img(:,:,2);
        end
        % Get binomial stats
        binom_st = CELLCOAV_BMRP(bw_vessel, ...
            cell_diam_um, umppix,tbl.tot_cells(ind),'n_obs_success',...
            tbl.obs_icf(ind)*tbl.tot_cells(ind));
        tbl.bmrp_icf(ind) = binom_st.binom_frac_mean;
        tbl.bmrp_cellcoav_p(ind) = binom_st.binom_p_nerandom;
        % Hypergeomtric model
        hyg_st = CELLCOAV_HMRP(bw_vessel, ...
            cell_diam_um, umppix,tbl.tot_cells(ind),'n_obs_success',...
            tbl.obs_icf(ind)*tbl.tot_cells(ind));
        tbl.hmrp_icf_mean(ind) = hyg_st.hmrp_icf_mean;
        tbl.hmrp_cellcoav_p(ind) = hyg_st.hmrp_cellcoav_p;
        % Fill in other metadata
        %             tbl.group_id(ind) = 2;
        tbl.sample_id(ind) = N;
        tbl.repl_id(ind) = n;
        
        
        writetable(tbl,[output_path '/arcas_results.csv']);
        waitbar(k/(grp_tot*N_tot*n_tot),hw);
        k=k+1;
    end
end

tbl.coloc_cells = tbl.obs_icf .* tbl.tot_cells;

close(hw);

% Export table
writetable(tbl,[output_path '/arcas_results.csv']);

% Read table
clear all
output_path = [getappdata(0,'proj_path') '/temp_data/StatsComparison'];
tbl = readtable([output_path '/arcas_results.csv']);
load([output_path '/arcas_results.mat']);


%Averaged Table
sample_tbl = grpstats(tbl,{'group_id','sample_id'},'mean','DataVars',...
    {'vld','tot_cells','coloc_cells','obs_icf','bmrp_icf','bmrp_cellcoav_pval'}) ;
grp1_ind = sample_tbl.group_id==1;
grp2_ind = sample_tbl.group_id==2;



% Group permutation testing
pval = ArcasGui_GroupPermutationTesting(tbl.hmrp_cellcoav_p, tbl.group_id,...
    output_path,1);
saveas(gcf,[output_path '/F_GroupPerm_comparison.fig']);

keyboard

fid=fopen([output_path '/Stats_Output.txt'],'w');

%% Pre-results stats
% Compare VLD
figure;
y=sample_tbl.mean_vld;
boxplot(y,sample_tbl.group_id);
set(gca,'fontsize',7,'fontname','helvetica');
set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
[h,p]=ttest2(y(grp1_ind),y(grp2_ind));
fprintf(fid,'A) VLD P Val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
    mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
    mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
    mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
ylabel('VLD (mm/mm2)');
yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
set(gcf,'position', [100 100 150 160]);
saveas(gcf,[output_path '/A_VLD_Comparison.fig']);

% Compare Total Cells
figure;
y=sample_tbl.mean_tot_cells;
boxplot(y,sample_tbl.group_id);
set(gca,'fontsize',7,'fontname','helvetica');
set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
[h,p]=ttest2(y(grp1_ind),y(grp2_ind));
fprintf(fid, 'B) Tot cells P Val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
    mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
    mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
    mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
ylabel('Cells/FOV');
yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
set(gcf,'position', [100 100 150 160]);
saveas(gcf,[output_path '/B_TotCells_Comparison.fig']);


%% Results stats
% Raw colocalization comparison
figure;
y=sample_tbl.mean_tot_cells .* sample_tbl.mean_obs_icf;
boxplot(y,sample_tbl.group_id);
set(gca,'fontsize',7,'fontname','helvetica')
set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
[h,p]=ttest2(y(grp1_ind),y(grp2_ind));
fprintf(fid, 'C) Coloc cells/FOV P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
    mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
    mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
    mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
ylabel('Coloc. Cells/FOV');
yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
set(gcf,'position', [100 100 150 160])
saveas(gcf,[output_path '/C_ColocCells-p-FOV_Comparison.fig']);



% coloc normalized to vld
figure;
y=sample_tbl.mean_coloc_cells./sample_tbl.mean_vld;
boxplot(y,sample_tbl.group_id);
set(gca,'fontsize',7,'fontname','helvetica')
set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
[h, p]=ttest2(y(grp1_ind),y(grp2_ind));
fprintf(fid, 'D) VLD Norm Coloc P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
    mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
    mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
    mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
ylabel('Coloc. Cells/Vess. Len. (mm)');
yl=ylim;xl=xlim;
axis([xl yl(1) yl(2)*1.2]);
set(gcf,'position', [100 100 150 160])
saveas(gcf,[output_path '/D_ColocCells-p-VL_Comparison.fig']);



% Colocalization normalized to injected cell number
figure;
y=sample_tbl.mean_obs_icf;
boxplot(y,sample_tbl.group_id);
set(gca,'fontsize',7,'fontname','helvetica')
set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
[h, p]=ttest2(y(grp1_ind),y(grp2_ind));
fprintf(fid, 'E) Fraction of Coloc cells P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
    mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
    mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
    mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
ylabel('Fraction of Coloc. Cells');
yl=ylim;xl=xlim;
axis([xl yl(1) yl(2)*1.2]);
set(gcf,'position', [100 100 150 160])
saveas(gcf,[output_path '/E_FracColocCells_Comparison.fig']);

fprintf(fid, 'F) GroupPerm Test: %f\n',pval);
fclose(fid);

save([output_path '/arcas_results.mat']);

keyboard
end

