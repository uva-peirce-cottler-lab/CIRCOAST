function [bw_kernel, rk, ck] = Generate_RestrictKernel(restrict_rad,vess_rad_pix)


% kernel image to restrict watershed spawn points
bw_kernel = false(3*restrict_rad);
ctr_pix = 1.5*(restrict_rad+1);
bw_kernel(ctr_pix,ctr_pix)=1;
kernel_ed = bwdist(bw_kernel);
bw_kernel = kernel_ed<restrict_rad;
% Horizontal Strip
bw_kernel(ctr_pix - ceil(4*vess_rad_pix):ctr_pix + ceil(2*vess_rad_pix),:)=1;
% Vertical Strip
bw_kernel(:, ctr_pix - ceil(4*vess_rad_pix):ctr_pix + ceil(2*vess_rad_pix))=1;
[rk, ck] =size(bw_kernel);



end