

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 14;
umppix = 424.5/1024;
vessel_rad_pix = 3;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = 20;
target_vld_mmpmm2=20;


imgs_path = 'D:\Box Sync\1. ARCAS\ExperimentData\Exp_259 ARCAS Imgs VEAP\veap';
imgs_items = dir([imgs_path '/*_bv.tif']);


for n=1:numel(imgs_items)
    
      bw = imread([imgs_path '/' imgs_items(n).name]);
      bw=bw(:,:,1);
      img =  imread([imgs_path '/' regexprep(imgs_items(n).name,'_bv','_proc')]);
      [mcm_colfrac_means(n), mcm_colfract_stds(n), ~, comp_img, ~] = ...
            ArcasGui_monteCarloSim_Driver(bw, cell_diam_um, umppix, ...
            tot_trials, tot_cells, 10000);
%         keyboard
       comp_img(:,:,2) = img(:,:,2);
       imwrite(comp_img, [imgs_path '/' regexprep(imgs_items(n).name,'_bv','_comp')]);
       
       
       imwrite(bwdist(bw)<cell_diam_um/umppix/2, ...
           [imgs_path '/' regexprep(imgs_items(n).name,'_bv','_celldil')]);
       
       comp_img(:,:,2) = im2uint8(bw).*.7;
       imwrite(comp_img, [imgs_path '/' regexprep(imgs_items(n).name,'_bv','_comp')]);
       
       bw_skel  = imread([imgs_path '/' regexprep(imgs_items(n).name,'_bv','_skel')]);
       vld(n) = vdl_pix2mmpmm2(bw_skel, umppix);
       vf(n) = sum(bw(:))./numel(bw);
end

f=fopen([imgs_path '/mcm_output.csv'],'w');
fprintf(f,'img_name,vf,vld_mmpmm2,mean_ccf,std_ccf,\n');
for n=1:numel(imgs_items)
   fprintf(f,'%s,%.6f,%.6f,%.6f,%.6f,\n', regexprep(imgs_items(n).name,'_bv',''),...
       vf(n),vld(n), mcm_colfrac_means(n),mcm_colfract_stds(n)); 
end
fclose(f);


