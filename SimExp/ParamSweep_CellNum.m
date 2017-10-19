clear all
close all

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 15;
umppix = 424.5/512;
vessel_rad_pix = 3;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = [1 10 25 50 100 500];


proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/CellNumber'];
if isempty(dir(out_path)); mkdir(out_path); end
delete([out_path '/*.tif']);

multiWaitbar('Images #', 0 );
multiWaitbar('Dilation #', 0 );
for r = 1:8
[bw_vessel, bw_skel_init, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, vessel_rad_pix,umppix,...
        restrict_rad, 'VesselLengthDensityLimits', [30 10]);
    vessel_ld_mmpmm2(r)=stat_st.VesselLengthDensity_mmpmm2;
    
    
     
    % Column: image dilations
    for c=1:numel(tot_cells)
        % Calculate mean col_frac
        [mcm_colfrac_means(r,c), mcm_colfract_stds(r,c), ~, comp_img, ~] = ...
            ArcasGui_monteCarloSim_Driver({bw_vessel}, cell_diam_um, umppix, ...
            tot_trials, tot_cells(c), 10000);
     
%             sbar_len = round(100/umppix);
%             comp_img(end-35:end-20,end-50-sbar_len:end-50,:)=...
%                 intmax(class(comp_img));
            comp_img(:,:,2) = uint8(comp_img(:,:,2)*.7);
            imwrite(comp_img, [out_path ...
                sprintf('/ParSweep_C%iR%i_cn%ium_VD%.4f_VLD%.4f.tif',...
                r,c,tot_cells(c),stat_st.VesselAreaFraction, stat_st.VesselLengthDensity_mmpmm2)]);
     
        
        % Calculate mean and std with binomial distribution
        binom_st = CELLCOAV_BMRP({bw_vessel}, cell_diam_um, umppix,tot_cells(c));
        bmd_colfrac_means(r,c) = binom_st.binom_frac_mean;
        bmd_colfrac_stds(r,c) = binom_st.binom_frac_std;
        
        
        
        
        % Randomly generate n binomial random fractions to match MCM
        brf = binornd(tot_cells(c),bmd_colfrac_means(r,c), [1 tot_trials])/tot_cells(c);
        bm_colfrac_means(r,c) = mean(brf);
        bm_colfrac_stds(r,c) = std(brf);
        
        %         keyboard
        
        
        multiWaitbar('Dilation #', c/6 );
    end
    multiWaitbar('Images #', r/8 );
end
multiWaitbar('Dilation #', 'Close' )
multiWaitbar('Images #', 'Close' )

save([out_path '/parsweep_cs_data.mat']);
load([out_path '/parsweep_cs_data.mat']);
% keyboard


x_data = 1:numel(tot_cells);
xtxt = 'Vessel Length Density (mm/mm^2)';
xspace =.0; 

figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),std(mcm_colfrac_means,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfrac_means,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar((1:numel(tot_cells)), mean(bm_colfrac_means,1),std(bm_colfrac_means,1),'rx');
% hE(3) = errorbar(x_data+xspace, mean(bmd_colfrac_means,1),std(bmd_colfrac_means,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_means,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
for n=1:numel(tot_cells)
       plot([x_data(n) x_data(n)],...
          [ya(1)*.90 mean(bm_colfrac_means(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. Mean')
xlabel('Cell Number')
% legend([hE(1) hE(2) hE(3)],{'MCM','BMRP-S','BMRP'})
 legend({'MCM'})
% x=xlim;y=ylim;
% axis([0 x(2) y(1) y(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
xticks(1:numel(tot_cells))
xticklabels(arrayfun(@(x) num2str(x),tot_cells,'UniformOutput',0));
axis([xa(1)-1 xa(2)+1 ya(1)*.90 ya(2)*1.1])
set(gcf,'Position', [100 100 260 260])
saveas(gcf,[out_path '/sweep_cn_mean_cell_coloc_fraction_mean.fig'])


figure('Units', 'pixels');
hold on; clear hE
hE(1) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),std(mcm_colfract_stds,1),'b.');
hE(2) = errorbar(x_data-xspace, mean(mcm_colfract_stds,1),zeros([1 numel(x_data)]),'b.');
% hE(2) = errorbar((1:numel(tot_cells)), mean(bm_colfrac_stds,1),std(bm_colfrac_stds,1),'rx');
% hE(3) = errorbar(x_data+xspace, mean(bmd_colfrac_stds,1),std(bmd_colfrac_stds,1),'r.');
% hE(4) = errorbar(x_data+xspace, mean(bm_colfrac_stds,1),zeros([1 numel(x_data)]),'r.');
xa=xlim;ya=ylim;
% axis([xa ya(1)*.90 ya(2)*1.1])
for n=1:numel(tot_cells)
   hold on
   plot([n n],[ya(1)*.90 mean(bm_colfrac_stds(:,n))],'Color',[.6 .6 .6],'LineStyle','--')
end
for n=1:numel(hE); hE(n).CapSize=6; end
hold off
ylabel('Mean of Cell Coloc. Fract. STD')
xlabel('Cell Number')
% legend({'MCM','BMRP-S','BMRP'})
% x=xlim;y=ylim;
% axis([0 x(2) y(1) y(2)])
beautifyAxis(gca);
set(gca, 'XGrid', 'off')
xticks(1:numel(tot_cells))
xticklabels(arrayfun(@(x) num2str(x),tot_cells,'UniformOutput',0));
axis([xa(1)-1 xa(2)+1 ya(1)*.90 ya(2)*1.1])
set(gcf,'Position', [100 100 260 260])
saveas(gcf,[out_path '/sweep_cn_mean_cell_coloc_fraction_std.fig'])



[RHO PVAL] = corr(tot_cells',mean(bm_colfrac_means,1)')