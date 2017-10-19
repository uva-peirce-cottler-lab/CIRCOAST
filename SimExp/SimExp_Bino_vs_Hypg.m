function [] = SimExp_Bino_vs_Hypg()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ClearPreviousData=1;
%Output Directory
output_path = [getappdata(0,'proj_path') '/temp_data/Bino Vs Hypg'];
mkdir(output_path);
delete([output_path '/*.*']);

% Simulation Parameters
N_tot=1000;
N_imgs = 50;
img_dim = [512 512];
umppix = 424.5/512;
cell_diam_um = 15;
cell_diam_pix = cell_diam_um/umppix;


% Data table 
% image_Name, group_id, sample_id, rep_id, vld, tot_cells, obs_col_cells, obs_icf, bmrp_icf
tbl=table;
tbl.img_name = cell(N_tot,1);
tbl.vld = zeros(N_tot,1);
tbl.bino_mean_icf = zeros(N_tot,1);
tbl.hypgs_mean_icf = zeros(N_tot,1);


% Requested vld is not exactly the same as actual
% Save time by regressing 4 additional times after each vessel network
% instead of producing from scratch
sub_vect = bsxfun(@times, rand(5,N_tot/5), [0; 5; 10; 15; 20]);
req_vld = normrnd(ones(N_imgs,1).*35,2);


% Create image dataset
hw = waitbar(0,'Processing Images');
bw_skel = false(img_dim);
k=1;
for N=1:N_imgs

               % Create Initial over dense vessel network
               [~, init_bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
                   3*6, 'VesselLengthDensityLimits', [50 0],'FigureVisible','on');
               % Identify central vessel to protect
               bw_protect = VesselGen_SelectProtectedSegment(init_bw_skel,'FigureVisible','on');

           
           % Iteratively remove vessel segments to a range of vessel densities
           [bw_vessel, bw_skel, stat_st] =...
                VesselGen_RegressNetwork(init_bw_skel, 3,umppix, ...
                'VesselLengthDensityTarget',req_vld(N),...
                'BW_Protected',bw_protect,'FigureVisible','on');
           
           
           for c = 1:N_tot/N_imgs
               tbl.vld(k) = stat_st.vld_mmpmm2;
               
               if N==1
                   % Distribute cells in image
                   [~, ~, ~,comp_img,tbl.whole_img_obs_icf(k)] = ...
                       CELLCOAV_MCMRP(bw_vessel, 1, ...
                       tbl.tot_cells(k), cell_diam_pix);
                   imwrite(comp_img,[output_path '/' tbl.img_name{N}]);
                   for n=1:4
                       imwrite(sub_img{n},[output_path '/' sprintf('I%i_%i.tif',N,n)]);
                   end
               end
               
                % Do a sweep of cell number values
                binom_st = CELLCOAV_BMRP(bw_vessel, ...
                cell_diam_um, umppix,tbl.tot_cells(N));
             	% Get p value form whole image
                tbl.bino_mean_icf(k) = binom_st.binom_frac_mean;
                tbl.bino_img_p_val(k) = binom_st.binom_p_nerandom;
                
                  % Do a sweep of cell number values
                binom_st = CELLCOAV_HMRP(bw_vessel, ...
                cell_diam_um, umppix,tbl.tot_cells(N));
             	% Get p value form whole image
                tbl.hypg_mean_icf(k) = binom_st.binom_frac_mean;
               
                keyboard
                k=k+1;
           end

           
         
            
            waitbar(k/(N_tot),hw);
            
 
end
 

close(hw);

% Export table
writetable(tbl,[output_path '/bino_vs_hypg_results.csv']);

% % Read table
% output_path = [getappdata(0,'proj_path') '/temp_data/StatsComparison'];
% tbl = readtable([output_path '/arcas_results.csv']);
% load([output_path '/arcas_results.mat']);
% 
% 
% %Averaged Table
% sample_tbl = grpstats(tbl,{'group_id','sample_id'},'mean','DataVars',...
%     {'vld','tot_cells','coloc_cells','obs_icf','bmrp_icf','bmrp_pval_gt_rand'}) ;
% grp1_ind = sample_tbl.group_id==1;
% grp2_ind = sample_tbl.group_id==2;
% 
% 
% 
% % Group permutation testing
% pval = ArcasGui_GroupPermutationTesting(tbl.bmrp_pval_gt_rand, tbl.group_id,...
%     output_path,1);
% saveas(gcf,[output_path '/F_GroupPerm_comparison.fig']);
% 
% keyboard
% 
% fid=fopen([output_path '/Stats_Output.txt'],'w');
% 
% %% Pre-results stats
% % Compare VLD
% figure;
% y=sample_tbl.mean_vld;
% boxplot(y,sample_tbl.group_id);
% set(gca,'fontsize',7,'fontname','helvetica');
% set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
% [h,p]=ttest2(y(grp1_ind),y(grp2_ind));
% fprintf(fid,'A) VLD P Val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
%     mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
%     mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
%     mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
% ylabel('VLD (mm/mm2)');
% yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
% set(gcf,'position', [100 100 150 160]);
% saveas(gcf,[output_path '/A_VLD_Comparison.fig']);
% 
% % Compare Total Cells
% figure;
% y=sample_tbl.mean_tot_cells;
% boxplot(y,sample_tbl.group_id);
% set(gca,'fontsize',7,'fontname','helvetica');
% set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
% [h,p]=ttest2(y(grp1_ind),y(grp2_ind));
% fprintf(fid, 'B) Tot cells P Val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
%     mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
%     mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
%     mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
% ylabel('Cells/FOV');
% yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
% set(gcf,'position', [100 100 150 160]);
% saveas(gcf,[output_path '/B_TotCells_Comparison.fig']);
% 
% 
% %% Results stats
% % Raw colocalization comparison
% figure;
% y=sample_tbl.mean_tot_cells .* sample_tbl.mean_obs_icf;
% boxplot(y,sample_tbl.group_id);
% set(gca,'fontsize',7,'fontname','helvetica')
% set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
% [h,p]=ttest2(y(grp1_ind),y(grp2_ind));
% fprintf(fid, 'C) Coloc cells/FOV P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
%     mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
%     mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
%     mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
% ylabel('Coloc. Cells/FOV');
% yl=ylim;xl=xlim; axis([xl yl(1) yl(2)*1.2]);
% set(gcf,'position', [100 100 150 160])
% saveas(gcf,[output_path '/C_ColocCells-p-FOV_Comparison.fig']);
% 
% 
% 
% % coloc normalized to vld
% figure;
% y=sample_tbl.mean_coloc_cells./sample_tbl.mean_vld;
% boxplot(y,sample_tbl.group_id);
% set(gca,'fontsize',7,'fontname','helvetica')
% set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
% [h, p]=ttest2(y(grp1_ind),y(grp2_ind));
% fprintf(fid, 'D) VLD Norm Coloc P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
%     mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
%     mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
%     mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
% ylabel('Coloc. Cells/Vess. Len. (mm)');
% yl=ylim;xl=xlim;
% axis([xl yl(1) yl(2)*1.2]);
% set(gcf,'position', [100 100 150 160])
% saveas(gcf,[output_path '/D_ColocCells-p-VL_Comparison.fig']);
% 
% 
% 
% % Colocalization normalized to injected cell number
% figure;
% y=sample_tbl.mean_obs_icf;
% boxplot(y,sample_tbl.group_id);
% set(gca,'fontsize',7,'fontname','helvetica')
% set(gca,'XTickLabel',{'Group1','Group2'},'fontsize',8);
% [h, p]=ttest2(y(grp1_ind),y(grp2_ind));
% fprintf(fid, 'E) Fraction of Coloc cells P val: %.2e, [%.2f+-%.2f, %.2f+-%.2f, %.1f]\n',p, ...
%     mean(y(sample_tbl.group_id==1)), std(y(sample_tbl.group_id==1)), ...
%     mean(y(sample_tbl.group_id==2)), std(y(sample_tbl.group_id==2)),...
%     mean(y(sample_tbl.group_id==2))./mean(y(sample_tbl.group_id==1))*100-100);
% ylabel('Fraction of Coloc. Cells');
% yl=ylim;xl=xlim;
% axis([xl yl(1) yl(2)*1.2]);
% set(gcf,'position', [100 100 150 160])
% saveas(gcf,[output_path '/E_FracColocCells_Comparison.fig']);
% 
% fprintf(fid, 'F) GroupPerm Test: %f\n',pval);
% fclose(fid);
% 
% save([output_path '/arcas_results.mat']);
% 
%     keyboard
% end
% 
