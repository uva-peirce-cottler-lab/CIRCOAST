clear all
close all

CELL_AREA_FRACT=0;
% dilCellFcn= @(x) ArcasGui_PerccentCellOverlap(imgs_cell,cell_diam_um,umppix)

% VF VLD u(MCM) std(MCM)
imgs_tbl = [0.229609 17.743761 0.485909 0.111664;...
            0.31171 24.456622 0.614923 0.108694;...
            0.110985 4.50438 0.197528 0.089027;...
            0.090047 6.8923 0.201033 0.089629;...
            0.525187 31.131524 0.797553 0.089943];
	

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 14;
umppix = 424.5/512;
vessel_rad_pix = 4;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = 20;
vessel_frac_target=.5;

proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/VLD_VF_WithReal'];
if isempty(dir(out_path)); mkdir(out_path); end

% Delete all files from previous analysis
delete([out_path '/*.tif']);



target_vessel_ld_mmpmm2 = [75 50 40 30 10 5];

fillFcn=@(bw, r) bwdist(bw) <=r;

multiWaitbar('Images #', 0 );
multiWaitbar('Dilation #', 0 );
for r = 1:8
    [bw_vessel, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, vessel_rad_pix,umppix,...
        restrict_rad, 'VesselLengthDensityLimits', [85 0]);
    
    % Remove all invisible hidden segments
    bw_skel = VesselGen_SimplifyNetwork(bw_skel);
    bw_skel = VesselGen_SimplifyNetwork(bw_skel);
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel);
    
    % Iteratively remove vessel segments to a range of vessel densities
    [~, skel_cell{r,1}, stat_st] =...
        VesselGen_RegressNetwork(bw_skel, vessel_rad_pix,umppix, ...
        'VesselLengthDensityTarget',target_vessel_ld_mmpmm2(1),'BW_Protected',bw_protect);
    obs_vessel_ld_mmpmm2(r,1) = stat_st.vd_mmpmm2;
    for c=2:6
        [~, skel_cell{r,c}, stat_st] =...
            VesselGen_RegressNetwork(skel_cell{r,c-1}, vessel_rad_pix,umppix, ...
            'VesselLengthDensityTarget',target_vessel_ld_mmpmm2(c),'BW_Protected',bw_protect);
        obs_vessel_ld_mmpmm2(r,c)=stat_st.vd_mmpmm2;
    end
    
    
    % Dilate each image fixed amount
    for c=1:6
        %Dilate image untilt he target vessel fraction is reached
        imgs_cell{r,c}=fillFcn(skel_cell{r,c}, vessel_rad_pix);
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
%            sbar_len = round(100/umppix);
%             comp_img(end-35:end-20,end-50-sbar_len:end-50,:)=...
%                 intmax(class(comp_img));
            
        comp_img(:,:,2) = uint8(comp_img(:,:,2)*.7);
        imwrite(comp_img, [out_path ...
            sprintf('/ParSweep_vf_C%iR%i_VD%.4f_VLD%.4f.tif',...
            r,c,vessel_frac_target, target_vessel_ld_mmpmm2(c))]);
         imwrite(im2uint8(cat(3,zeros(img_dim),comp_img(:,:,2),zeros(img_dim))), [out_path ...
            sprintf('/ParSweep_vf_C%iR%i_VD%.4f_VLD%.4f_vessel.tif',...
            r,c,vessel_frac_target, target_vessel_ld_mmpmm2(c))]);
        
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

xtxt = 'Vessel Fraction';
x_data = mean(obs_vessel_frac,1);
xspace=0;


keyboard






% range_x_data = max(obs_vessel_frac,[],1)-  min(obs_vessel_frac,[],1);
figure('Units', 'pixels');
hold on
hE(1) = errorbar(x_data, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),zeros([1 numel(x_data)]),'b.');
hP(1) = plot(imgs_tbl(:,1),imgs_tbl(:,3),'rx');
yr=[0.02 0.04];
for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [0 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean    ')
xlabel(xtxt)
% legend([hE(1) hP(1)],{'MCM','Exp'})
xa=xlim;ya=ylim;
axis([xa(1) xa(2) ya(1) ya(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 500 300])
saveas(gcf,[out_path '/sweep_vf_mean_cell_coloc_fraction_mean.fig'])




% Plotting
x_data = target_vessel_ld_mmpmm2;
xtxt = 'Vessel Length Density (mm/mm^2)';
figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),zeros([1 numel(x_data)]),'b.');
hP(1) = plot(imgs_tbl(:,2),imgs_tbl(:,3),'rx');
for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [0 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean   ')
xlabel(xtxt)
% legend([hE(1) hP(1)],{'MCM','Exp'})
xa=xlim;ya=ylim;
axis([xa(1) xa(2) ya(1) ya(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 500 300])
saveas(gcf,[out_path '/sweep_vld_mean_cell_coloc_fraction_mean.fig'])
