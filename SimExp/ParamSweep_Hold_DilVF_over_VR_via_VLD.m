

clear all
close all

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 15;
umppix = 424.5/512;
vessel_rad_um = [2 4 8 12 15 20];
vessel_rad_pix = vessel_rad_um.*umppix;
restrict_rad = 3*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = 20;

cell_dil_vf_target = .40;




proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/Hold_DilVF_over_VR_via_VLD'];
if isempty(dir(out_path)); mkdir(out_path); end
delete([out_path '/*.*']);

multiWaitbar('Images #', 0 );
multiWaitbar('Dilation #', 0 );
for r = 1:8
    [bw_vessel, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
        restrict_rad, 'VesselLengthDensityLimit', [40 0]);
    %     vessel_ld_mmpmm2(r)=stat_st.VesselLengthDensity_mmpmm2;
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel);
    
    bw_skel2=VesselGen_RegressNetwork(bw_skel,3,umppix, ...
        'VesselLengthDensityTarget',25,'BW_Protected',bw_protect);
    
    
    % Regress skeleton for each vessel img until DilVF is constant
    for c=1:6
        [imgs_cell{r,c}, ~, stat_st] =...
            VesselGen_RegressNetwork(bw_skel2,vessel_rad_pix(c),umppix, ...
            'Cell_Rad_Pix', cell_diam_um/umppix/2,...
            'VesselAreaFractionTarget',cell_dil_vf_target ,'BW_Protected',bw_protect);
        obs_vld_mmpmm2(r,c)=stat_st.vd_mmpmm2;
        obs_vf(r,c)=stat_st.vf;
        obs_binom_vf(r,c) = CELLCOAV_CellDilateVessFrac(imgs_cell(r,c),cell_diam_um,umppix);
        
        obs_dil_vf(r,c) = sum(sum(bwdist(imgs_cell{r,c}) <= ...
            cell_diam_um/umppix/2))./numel(imgs_cell{r,c});
    end
    %     keyboard
    %      for c=1:6
    %     subplot(2,3,c)
    %     imshow(imgs_cell{r,c})
    %      end
    multiWaitbar('Dilation #', 0 );
    
    % Column: image dilations
    for c=1:6
        % Calculate mean col_frac
        [mcm_colfrac_means(r,c), mcm_colfract_stds(r,c), ~, comp_img, ~] = ...
            ArcasGui_monteCarloSim_Driver(imgs_cell(r,c), cell_diam_um, umppix, ...
            tot_trials, tot_cells, 10000);
        
        bw_dil_area = bwdist(comp_img(:,:,2))<=cell_diam_um & ~comp_img(:,:,2);
        bw_dil_area= bw_dil_area + imbinarize(comp_img(:,:,3));
        %         sbar_len = round(100/umppix);
        comp_img(:,:,3)=im2uint8(bw_dil_area);
        %         comp_img(end-35:end-20,end-50-sbar_len:end-50,:)=...
        %             intmax(class(comp_img));
        comp_img(:,:,2) = uint8(comp_img(:,:,2)*.7);
        imwrite(comp_img, [out_path ...
            sprintf('/ParSweep_vf_C%iR%i_VR%.4f_VLD%.4f.tif',...
            r,c,vessel_rad_um(c), obs_vld_mmpmm2(r,c))]);
        
        
        % Calculate mean and std with binomial distribution
        binom_st = CELLCOAV_BMRP(imgs_cell(r,c), cell_diam_um, umppix,tot_cells);
        bmd_colfrac_means(r,c) = binom_st.binom_frac_mean;
        bmd_colfrac_stds(r,c) = binom_st.binom_frac_std;
        
        
        
        
        % Randomly generate n binomial random fractions to match MCM
        brf = binornd(tot_cells,bmd_colfrac_means(r,c), [1 tot_trials])/tot_cells;
        bm_colfrac_means(r,c) = mean(brf);
        bm_colfrac_stds(r,c) = std(brf);
        
        %         keyboard
        
        
        multiWaitbar('Dilation #', c/6 );
    end
    multiWaitbar('Images #', r/8 );
end
multiWaitbar('Dilation #', 'Close' )
multiWaitbar('Images #', 'Close' )

save([out_path '/parsweep_data.mat']);
load([out_path '/parsweep_data.mat']);
keyboard

x_data = vessel_rad_um;
xtxt = 'Vessel Radius (um)';
xspace =0.0;

[p,tbl,stats] = anova1(mcm_colfrac_means);


figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar((1:numel(tot_cells)), mean(bm_colfrac_means,1),std(bm_colfrac_means,1),'rx');
%  hE(1) = errorbar(x_data+xspace, mean(bmd_colfrac_means,1),std(bmd_colfrac_means,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_means,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
axis([xa ya(1)*.90 ya(2)*1.1])
for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [ya(1)*.90 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean        ')
xlabel('Vessel Radius (um)')
% legend([hE(1) hE(3)],{'MCM','BMRP'})
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 210 210])
saveas(gcf,[out_path '/sweep_mean_cell_coloc_fraction_mean.fig'])



figure('Units', 'pixels');
hold on
hE(1) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),std(mcm_colfract_stds,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar((1:numel(tot_cells)), mean(bm_colfrac_stds,1),std(bm_colfrac_stds,1),'rx');
% hE(3) = errorbar(x_data+xspace, mean(bmd_colfrac_stds,1),std(bmd_colfrac_stds,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_stds,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
axis([xa ya(1)*.98 ya(2)*1.02])
for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [ya(1)*.98 mean(bm_colfrac_stds(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. STD')
xlabel('Vessel Radius (um)')
% legend({'MCM','BMRP-S','BMRP'})
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 240 240])
saveas(gcf,[out_path '/sweep_mean_cell_coloc_fraction_std.fig'])


[RHO,PVAL] = corr(x_data',mean(bmd_colfrac_means,1)')

