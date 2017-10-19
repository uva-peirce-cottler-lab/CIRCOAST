function [p_success, npix_tot] = CELLCOAV_CellDilateVessFrac(imgs_cell,cell_diam_um,umppix)
if ~iscell(imgs_cell); imgs_cell = {imgs_cell}; end

for n = 1:numel(imgs_cell)
    %         cell_diam_um/umppix
    %         bw_dil = imdilate(imgs_cell{n},strel('disk',ceil((cell_diam_um/umppix)/2)));
    
    % Each pixel distance to closest forground
    ed = bwdist(imgs_cell{n});
    
    % Calculate minimum distance from each cell to foreground
    bw_dil = ~((ed - (cell_diam_um/umppix)/2)>0);
    %         keyboard
    npix_true(n) = sum(bw_dil(:));
    npix(n) = numel(bw_dil);
end
npix_tot=sum(npix);
p_success = sum(npix_true)/npix_tot;

