clear all
close all

% Calculate mean coloc fraction over range of cecl


% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 15;
umppix = 424.5/512;
vessel_rad_pix = 3;
restrict_rad = vessel_rad_pix*6;
img_dim = [512 512];
tot_trials = 1e6;
tot_cells = 20;
target_vld_mmpmm2=20;

% vessel_frac_target=.65;

proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/Binom_Illustration'];
if isempty(dir(out_path)); mkdir(out_path); end
delete([out_path '/*.*']);



    [bw_vessel_init, bw_skel_init, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, vessel_rad_pix,umppix,...
        restrict_rad, 'VesselLengthDensityLimits', [35 10]);
    
    sbar_len = round(25/umppix)*512/114;
    
    imtool(bw_vessel_init)
    
    img1=zeros([512 512 3],'uint8');
    img1(:,:,2) = imresize(im2uint8(BW),[512 512],'nearest');
    img1(end-35:end-20,end-50-sbar_len:end-50,:)=...
                intmax(class(img1));
    figure;imshow(img1)
    
    img2=zeros([512 512 3],'uint8');
    rs_bw = imresize(im2uint8(BW),[512 512],'nearest');
    img2(:,:,2) = rs_bw;
    cell_dill_gs = im2uint8(bwdist(BW) <= (cell_diam_um/umppix/2) & ~BW);
    img2(:,:,3)=imresize(cell_dill_gs,[512 512], 'nearest');
      img2(end-35:end-20,end-50-sbar_len:end-50,:)=...
                intmax(class(img2));
    figure;imshow(img2)
    
  
    onv = [34 36; 31 69;54 87];
    offv = [95 88;68 26];
    bw_on= false(size(BW));
    bw_on(sub2ind(size(BW),onv(:,1),onv(:,2)))=1;
    bw_on = bwdist(bw_on) <=(cell_diam_um/umppix/2-.3);
    
    bw_off = false(size(BW));
    bw_off(sub2ind(size(BW),offv(:,1),offv(:,2)))=1;
    bw_off = bwdist(bw_off) <=(cell_diam_um/umppix/2-.3);
    
    on_gs = im2uint8(imresize(bw_on,[512 512],'nearest'));
    off_gs = im2uint8(imresize(bw_off,[512 512],'nearest'));
    ;
    img3 = img2;
    img3(:,:,1)=on_gs;
    temp=img3(:,:,1);
    temp(off_gs>0)=155;
    img3(:,:,1)=temp;
    img3(end-35:end-20,end-50-sbar_len:end-50,:)=...
                intmax(class(img3));
    
    imwrite(img1,[out_path '/img1.tif'])
    imwrite(img2,[out_path '/img2.tif'])
    imwrite(img3,[out_path '/img3.tif'])
    
    save([out_path '/data.mat'])
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel_init);
    
    % Generate Skeleton Image
    [bw_vessel{r,1},~, stat_st] =...
        VesselGen_RegressNetwork(bw_skel_init, vessel_rad_pix,umppix, ...
        'VesselLengthDensityTarget',target_vld_mmpmm2,'BW_Protected',bw_protect);
    obs_vessel_ld_mmpmm2(r) = stat_st.vd_mmpmm2;
    
    