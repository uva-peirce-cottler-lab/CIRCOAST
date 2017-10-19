function vld_mmpmm2= ImageQuant_Pix2VesselLengthDensity(bw_skel,umppix)
img_dim = size(bw_skel);

vld_mmpmm2 =  sum(bw_skel(:))*(umppix/1000)/(img_dim(1)*umppix/1000).^2;