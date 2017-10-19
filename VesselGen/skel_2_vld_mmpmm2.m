function vld = skel_2_vld_mmpmm2(bw_skel, umppix)
img_dim = size(bw_skel);

vld = sum(bw_skel(:))*(umppix/1000)/(img_dim(1)*umppix/1000).^2;