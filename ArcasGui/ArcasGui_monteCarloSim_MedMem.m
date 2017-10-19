function [hit_mean, hit_std, cell_vasc_dists, comp_image, frac_hits, ...
    frac_cell_densities, frac_cell_overlaps] = ...
    ArcasGui_monteCarloSim_MedMem(img_in, num_trials, num_events, event_diam)
% Reads in a black and white mask image ('image_IN') to restrict where
% coordinates can be placed. In other words, 'image_IN" needs to be binary.

%'num_events' is the number of cells in the simulation

% 'pixelsize' is the size of cells (diameter in pixels)

% 'AVG' is the percent coverage of cells in the randomized simulation.
% tic

assert(mod(event_diam,1)==0, 'event diam must to be an integer');

if ischar(img_in)
img = im2double((imread(img_in)));
else img=img_in;
end

[r, c] = size(img);

% Tool used to make each cell center into a blob
cell_strel = strel('disk', ceil((event_diam+1)/2));

% Pericyte center matrix
cell_matrix = false(r,c);

% Randomly create cell centers
cell_center_ind = randi([1 r*c],num_trials, num_events);

% Offset each row by 1024^2 to shift to next z slice
% offset = cumsum(1024^2*ones(num_trials,1))-1024^2;
offset = cumsum(numel(img)*ones(num_trials,1))-numel(img);

% Convert indices for peri center so each row of peri_ctr_ind is for a
% seperate z slice of the peri_zimg
cell_ctr_zstack_ind = bsxfun(@plus,offset,cell_center_ind);

% Z stack of center of each cell
cell_ctr_zimg = false([size(img) num_trials]);
cell_ctr_zimg(cell_ctr_zstack_ind(:)) = 1;

% For each of the cell images, expand center point to 2D blob
cell_zimg = imdilate(cell_ctr_zimg, cell_strel);

% Euclidian distance transform input img
ed = bwdist(img);

cell_vasc_dists_cell = cell(1,num_trials);
for n=1:num_trials

    % Connected component analysis - looks for continuous blobs and records
    % the list of pixel indeces that make up that blob
    cc = bwconncomp(cell_zimg(:,:,n));
    
    % Go through each blob, grab distance values, find min
    for k = 1:cc.NumObjects
        dist_values = ed(cc.PixelIdxList{k});
        cell_min_dist2_vasculature(k) = min(dist_values);
    end
    
    % Calc number of cells that hit and %
    current_hit_number = sum(cell_min_dist2_vasculature==0);
    current_percent_hit = current_hit_number/numel(cell_min_dist2_vasculature);
    
    frac_hits(n,1)= current_percent_hit; %this is where i put in the percent coverage from one loop
    
    % Fraction of Image area that is a cell
    frac_cell_densities(n,1) = sum(sum(cell_zimg(:,:,n)))/(r*c);
    
    % Fraction of cells that overlap
    frac_cell_overlaps(n,1) = (num_events-cc.NumObjects)/num_events;
    
    % Clear 
    cell_matrix(:) = 0;
    
    cell_vasc_dists_cell{n} = cell_min_dist2_vasculature;
end

% For each perictye across all trials, min distance to vasculature
cell_vasc_dists = horzcat(cell_vasc_dists_cell{:});

% Make composite image with cells red
comp_image = im2uint8(cat(3, img,img,img));
temp = comp_image(:,:,1); temp(cell_zimg(:,:,end)) = 256; comp_image(:,:,1) = temp;
temp = comp_image(:,:,2); temp(cell_zimg(:,:,end)) = 0;   comp_image(:,:,2) = temp;
temp = comp_image(:,:,3); temp(cell_zimg(:,:,end)) = 0;   comp_image(:,:,3) = temp;

% Return stats for the trials run
hit_mean = mean(frac_hits)*100;
hit_std = std(frac_hits);

% figure; imshow(comp_image);
keyboard 
