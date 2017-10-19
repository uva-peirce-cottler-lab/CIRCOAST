clear all
close all

% Calculate mean coloc fraction over range of cecl


% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = [5 10 15 20 30 40];
umppix = 424.5/512;
vessel_rad_pix = 3;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = 20;
target_vld_mmpmm2=20;

% vessel_frac_target=.65;

proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/CellDiam'];
if isempty(dir(out_path)); mkdir(out_path); end
delete([out_path '/*.*']);


% fillFcn=@(bw, r) imdilate(bw,strel('disk',r,0)) | imdilate(bwmorph(bw,'branchpoints'),strel('disk',ceil(r*1.25),0));

multiWaitbar('Images #', 0 );
multiWaitbar('Dilation #', 0 );
for r = 1:8
    [bw_vessel_init, bw_skel_init, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, vessel_rad_pix,umppix,...
        restrict_rad, 'VesselLengthDensityLimits', [30 0]);
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel_init);
    
    % Generate Skeleton Image
    [bw_vessel{r,1},~, stat_st] =...
        VesselGen_RegressNetwork(bw_skel_init, vessel_rad_pix,umppix, ...
        'VesselLengthDensityTarget',target_vld_mmpmm2,'BW_Protected',bw_protect);
    obs_vessel_ld_mmpmm2(r) = stat_st.vd_mmpmm2;
    
    
    multiWaitbar('Dilation #', 0 );
    % Calculate Metrics from images
    % Column: image dilations
    for c=1:6
        % Calculate mean col_frac
        [mcm_colfrac_means(r,c), mcm_colfract_stds(r,c), ~, comp_img, ~] = ...
            ArcasGui_monteCarloSim_Driver(bw_vessel{r,1}, cell_diam_um(c), umppix, ...
            tot_trials, tot_cells, 10000);
%         bw_dil_area = bwdist(comp_img(:,:,2))<=cell_diam_um(c) & ~comp_img(:,:,2);
%         bw_dil_area= bw_dil_area + imbinarize(comp_img(:,:,3));
%         sbar_len = round(100/umppix);
%         comp_img(:,:,3)=im2uint8(bw_dil_area);
%         comp_img(end-35:end-20,end-50-sbar_len:end-50,:)=...
%             intmax(class(comp_img));
        comp_img(:,:,2) = uint8(comp_img(:,:,2)*.7);
        imwrite(comp_img, [out_path ...
            sprintf('/ParSweep_vf_C%iR%i_CS%.4f_VLD%.4f.tif',...
            r,c,cell_diam_um(c), target_vld_mmpmm2)]);
        
        
        % Calculate mean and std with binomial distribution
        binom_st = CELLCOAV_BMRP(bw_vessel{r,1}, cell_diam_um(c), umppix,tot_cells);
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
% keyboard
  
 
x_data = cell_diam_um;
xtxt = 'Cell Diameter (um)';
xspace = 0;


figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar(x_data+.5, mean(bmd_colfrac_means,1),std(bmd_colfrac_means,1),'ko');
% hE(3) = errorbar(x_data+xspace, mean(bm_colfrac_means,1),std(bm_colfrac_means,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_means,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
axis([x_data(1)-2 x_data(end)+2 ya(1)*.99 ya(2)*1.01])
for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [ya(1)*.99 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean')
xlabel(xtxt)
% legend([hE(1) hE(2) hE(3)],{'MCM','BMRP'})
% legend([hE(1) hE(3)],{'MCM','BMRP'})
beautifyAxis(gca);
set(gca,'TickLength'  , [0 .02])
set(gca, 'XGrid', 'off')
set(gcf,'Position', [100 100 260 260])
saveas(gcf,[out_path '/sweep_vld_mean_cell_coloc_fraction_mean.fig'])


figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),std(mcm_colfract_stds,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar(x_data+.5, mean(bmd_colfrac_stds,1),std(bmd_colfrac_stds,1),'ko');
% hE(3) = errorbar(x_data+xspace, mean(bm_colfrac_stds,1),std(bm_colfrac_stds,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_stds,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
axis([x_data(1)-2 x_data(end)+2 ya(1)*.99 ya(2)*1.01])

for n=1:numel(x_data)
    plot([x_data(n) x_data(n)],...
        [ya(1)*.95 mean(bm_colfrac_stds(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. STD')
xlabel(xtxt)
% legend([hE(1) hE(3)],{'MCM','BMRP'})
xa=xlim;ya=ylim;
axis([xa(1) xa(2) ya(1) ya(2)])
beautifyAxis(gca);
set(gca,'TickLength'  , [0 .02])
set(gcf,'Position', [100 100 260 260])
saveas(gcf,[out_path '/sweep_vf_mean_cell_coloc_fraction_std.fig'])

[RHO PVAL] = corr(x_data',mean(bm_colfrac_means,1)')

