function [hit_mean, hit_std, cell_dists_vect, comp_image, frac_hits, var_mem] = ...
    CELLCOAV_MCMRP(imgs_cell, num_trials, num_events, event_diam)
% Reads in a black and white mask image ('image_IN') to restrict where
% coordinates can be placed. In other words, 'image_IN" needs to be binary.
%'num_events' is the number of cells in the simulation
% 'pixelsize' is the size of cells (diameter in pixels)
% 'AVG' is the percent coverage of cells in the randomized simulation.
% tic

if (isnumeric(imgs_cell) || islogical(imgs_cell)) && all(size(imgs_cell>1));
imgs_cell = {imgs_cell};    
end
 
% assert(mod(event_diam,1)==0, 'cell diameter must to be an integer');

[bw_imgs, bw_rois] = ArcasGui_CreateHorizontal_ConcatImage(imgs_cell);

% Indiceies of valid spawn points
ind = reshape(1:numel(bw_rois),size(bw_rois));
spawn_indicies = ind(bw_rois);


% Randomly create cell centers
random_ints = randi([1 numel(spawn_indicies)], [num_trials num_events]);
cell_center_ind = spawn_indicies(random_ints);
if num_trials==1 && size(cell_center_ind,1)~=1; cell_center_ind=cell_center_ind'; end
 
% Test that indicies are within range of images
% bw_test = false(size(bw_imgs));
% bw_test(cell_center_ind)=1;
% imshow(bw_test)

% Each pixel distance to closest forground
ed = bwdist(bw_imgs);

 
% Calculate minimum distance from each cell to foreground
cell_dists = ed(cell_center_ind) - event_diam/2;
frac_hits = sum(cell_dists <=0,2)/num_events;

% Output arg vectorized results
cell_dists_vect = cell_dists(:)';


% Make composite image with cells red
% Green: Blood vessel
% Bright Red: On vessel Cell
% Dark Red: off Vessel Cell
bw_vessel=bw_imgs;
bw_cell_ctr = zeros(size(bw_vessel),'uint8');
bw_cell_ctr(cell_center_ind(1,:)) = 1+(cell_dists(1,:)<0);
% bw_cell_ctr(cell_center_ind(1,:))=1;
bw_cell = bwdist(bw_cell_ctr)<=event_diam/2;
% keyboard
bw_cell_ov = (bw_cell & bw_vessel) + bw_cell;
cc_ov = bwconncomp(bw_cell_ov);
gs_cell =zeros(size(bw_imgs),'uint8');
for n=1:cc_ov.NumObjects
%     fprintf('IsOV: %i\n',any(bw_cell_ov(cc_ov.PixelIdxList{n})==2))
    gs_cell(cc_ov.PixelIdxList{n})= any(bw_cell_ov(cc_ov.PixelIdxList{n})==2)*70+170;
end


comp_image = cat(3,im2uint8(gs_cell),im2uint8(bw_vessel),bw_cell_ctr);


% Return stats for the trials run
hit_mean = mean(frac_hits)*100;
hit_std = std(frac_hits);


var_mem = sum(arrayfun(@(x) x.bytes, whos));

% keyboard
