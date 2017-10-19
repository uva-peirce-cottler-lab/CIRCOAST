clear all
close all

CELL_AREA_FRACT=1;
% dilCellFcn= @(x) ArcasGui_PerccentCellOverlap(imgs_cell,cell_diam_um,umppix)

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 10; 
umppix = 424.5/512;
vessel_rad_pix = 3;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6; 
tot_cells = 20;
vessel_frac_target=.5;

proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/VF_via_VLD_Hold_DilVF'];
if isempty(dir(out_path)); mkdir(out_path); end
delete([out_path '/*.tif']);


% Make set of images with same vessel radiu
% Make a set of images with range of VF by
% Adjust VLD so that CelVF is held constant


fillFcn=@(bw, r) imdilate(bw,strel('disk',r,0)) | imdilate(bwmorph(bw,'branchpoints'),strel('disk',ceil(r*1.25),0));
   
multiWaitbar('Images #', 0 );
multiWaitbar('Dilation #', 0 );
for r = 1:8 
    [bw_vessel, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, vessel_rad_pix,umppix,...
        restrict_rad, 'VesselLengthDensityLimits', [55 0]);
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel);
 
    % Iteratively remove vessel segments to a range of vessel densities
     [~, bw_rmv_skel , stat_st] =...
            VesselGen_RegressNetwork(bw_skel, vessel_rad_pix,umppix, ...
            'VesselLengthDensityTarget',45,'BW_Protected',bw_protect);
        
      for c=1:6 
      vessel_cell{r,c} = bwdist(bw_rmv_skel)< x*2
    end
    
        
     obs_vessel_ld_mmpmm2(r,1) = stat_st.vd_mmpmm2;\
     obs_vf(r,1) = stat_st.vf;
    for c=2:6
        [~, skel_cell{r,c}, stat_st] =...
            VesselGen_RegressNetwork(skel_cell{r,c-1}, vessel_rad_pix,umppix, ...
            'VesselLengthDensityTarget',target_vessel_ld_mmpmm2(c),'BW_Protected',bw_protect);
        obs_vessel_ld_mmpmm2(r,c)=stat_st.vd_mmpmm2;
        
    end
    
    
    % Dilate each image to reach same vessel frac value
    for c=1:6 
        %Dilate image untilt he target vessel fraction is reached
     
        [imgs_cell{r,c} vessel_frac(r,c)]= ...
            VesselGen_Dilate2VesselFractionTarget(skel_cell{r,c}, vessel_frac_target,...
            cell_diam_um,umppix,CELL_AREA_FRACT);
        obs_vessel_frac(r,c) = sum(sum(imgs_cell{r,c}))./numel(imgs_cell{r,c}); 
    end
    
    figure; for c=1:6; subplot(2,3,c);imshow(imgs_cell{r,c}); end
    pause(2); close(gcf)
    
    multiWaitbar('Dilation #', 0 );
    % Calculate Metrics from images
    % Column: image dilations
    for c=1:6
        % Calculate mean col_frac
        [mcm_colfrac_means(r,c), mcm_colfract_stds(r,c), ~, comp_img, ~] = ...
            ArcasGui_monteCarloSim_Driver(imgs_cell(r,c), cell_diam_um, umppix, ...
            tot_trials, tot_cells, 10000);
        if r>0
            sbar_len = round(100/umppix);
            comp_img(end-35:end-20,end-50-sbar_len:end-50,:)=...
                intmax(class(comp_img));
        
            imwrite(comp_img, [out_path ...
                sprintf('/ParSweep_vf_C%iR%i_VD%.4f_VLD%.4f.tif',...
                r,c,vessel_frac_target, target_vessel_ld_mmpmm2(c))]);
        end
        
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

save([out_path '/parsweep_vld_data.mat']);
load([out_path '/parsweep_vld_data.mat']);
keyboard


x_data = target_vessel_ld_mmpmm2;
xtxt = 'Vessel Length Density (mm/mm^2)';

figure('Units', 'pixels');
hold on
hE(1) = errorbar(x_data-2, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'bx');
hE(2) = errorbar(x_data, mean(bm_colfrac_means,1),std(bm_colfrac_means,1),'rx');
hE(3) = errorbar(x_data+2, mean(bmd_colfrac_means,1),std(bmd_colfrac_means,1),'ko');
for n=1:numel(x_data)
      plot([x_data(n) x_data(n)],...
          [0 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=8; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean')
xlabel(xtxt)
legend([hE(1) hE(2) hE(3)],{'MCM','BMRP-S','BMRP'})
xa=xlim;ya=ylim;
axis([xa(1) xa(2) ya(1) ya(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 280 280])
saveas(gcf,[out_path '/sweep_vld_mean_cell_coloc_fraction_mean.fig'])


figure('Units', 'pixels');
hold on
hE(1) = errorbar(x_data-2, mean(mcm_colfract_stds,1),std(mcm_colfract_stds,1),'bx');
hE(2) = errorbar(x_data, mean(bm_colfrac_stds,1),std(bm_colfrac_stds,1),'rx');
hE(3) = errorbar(x_data+2, mean(bmd_colfrac_stds,1),std(bmd_colfrac_stds,1),'ko');
for n=1:numel(x_data)
      plot([x_data(n) x_data(n)],...
          [0 mean(bm_colfrac_stds(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=8; end
hold off
ylabel('Mean of Cell Coloc. Fract. STD')
xlabel(xtxt)
legend({'MCM','BMRP-S','BMRP'})
xa=xlim;ya=ylim;
axis([xa(1) xa(2) ya(1) ya(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 280 280])
saveas(gcf,[out_path '/sweep_vf_mean_cell_coloc_fraction_std.fig'])
