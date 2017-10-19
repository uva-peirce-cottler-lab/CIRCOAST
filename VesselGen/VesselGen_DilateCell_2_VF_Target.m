function [cell_rad_pix vf] = ...
    VesselGen_DilateCell_2_VF_Target(bw_vessel, target_vf)
DBUG=0;


cell_rad_pixs(1) = 1;
for n=1:1e8
    
   cell_dil_vf(n) = sum(sum(bwdist(bw_vessel) <= cell_rad_pixs(n)))./numel(bw_vessel);
   
   if cell_dil_vf(n)>= target_vf;
%       keyboard 
         [~,ind] = min([abs(cell_dil_vf(n)-target_vf) abs(cell_dil_vf(n-1)-target_vf)]);
         best_n = n-(ind-1);
         cell_rad_pix = cell_rad_pixs(best_n);
         vf = cell_dil_vf(best_n);
         break
   end
    
    cell_rad_pixs(n+1)= cell_rad_pixs(n)+0.1;
end