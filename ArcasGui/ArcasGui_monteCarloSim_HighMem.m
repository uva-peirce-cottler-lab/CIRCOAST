function [hit_mean, hit_std, comp_image] = ArcasGui_monteCarloSim_HighMem(img_in,...
    num_trials, num_events)
% Reads in a black and white mask image ('img_in') to restrict where
% coordinates can be placed. In other words, 'img_in" needs to be binary.

%'num_events' is the number of cells in the simulation

% 'pixelsize' is the size of cells (diameter in pixels)

% 'AVG' is the percent coverage of cells in the randomized simulation.
tic
% profile on
% Import mask image (best to use .tif)

if ischar(img_in)
J = im2double((imread(img_in)));
else J=img_in;
end

% q=zeros(10,1);
[r, c] = size(J);

% Tool used to make each cell center into a blob
cell_strel = strel('disk', 9);
cell_cc_strel = zeros(3,3,3);
cell_cc_strel(:,:,2)=1;

% Pericyte center matrix
peri_matrix = false(r,c);

% Randomly create cell centers
peri_center_ind = randi([1 r*c],num_trials, num_events);

% Offset each row by 1024^2 to shift to seperate z slices
offset = cumsum(1024^2*ones(num_trials,1))-1024^2;

% Convert indices for peri center so each row of peri_ctr_ind is for a
% seperate z slice of the peri_zimg
peri_ctr_zstack_ind = bsxfun(@plus,offset,peri_center_ind);

% Pericyte Z stack
peri_zimg = false([size(J) num_trials]);
peri_zimg(peri_ctr_zstack_ind(:)) = 1;
% squeeze(sum(sum(peri_zimg,1),2))


% For each of the cell images, expand center point to 2D blob
cell_zimg = imdilate(peri_zimg, cell_strel);

% Euclidian distance transform input img
ed = bwdist(J);

% Seperate blobs on each zslices indepedently (no blobs across zslices)
cc = bwconncomp(cell_zimg, cell_cc_strel);

% Test if connectivity kept each slice seperate, did not merge blobs across
% z = unique(cellfun(@(x) numel(unique(floor(x/1024^2+1))), cc.PixelIdxList));

% Calc which z slice each cell is from
z = cellfun(@(x) unique(floor(x/1024^2+1)), cc.PixelIdxList)';
z_matrix = bsxfun(@plus, z, zeros(size(z,1),num_trials));
unq_z_matrix = bsxfun(@eq, 1:num_trials, z_matrix); % 0 past 10000

% Obtain min distance from each blob
min_dist = cellfun(@(x) min(ed(mod(x,r*c)+1)),cc.PixelIdxList)';
min_dist_matrix = bsxfun(@plus, min_dist, zeros(size(min_dist,1),num_trials));

% Element >1 is min distance to vasculature
% Element of 1 means a hit
% Element of zero means that blob belonged to another z slice
unq_min_dist_matrix = (min_dist_matrix+1).*unq_z_matrix;

% Hit number for each trial is # of blobs with zero distance (1 in this case
% because we offset it)
hit_number = sum(unq_min_dist_matrix == 1);
% Number of cells in each trial (zslice), sometimes cells can
% overlap, be counted as one
cell_number = sum(unq_min_dist_matrix~=0,1);
% Average min distance of each cell from vasulcature
average_min_distance = sum(unq_min_dist_matrix,1)./cell_number;
% Percentage of cells that hit for each trial
cell_hit_prc = hit_number./cell_number;


% Make composite image with cells red
comp_image = cat(3, J,J,J);
temp = comp_image(:,:,1); temp(cell_zimg(:,:,end)) = 1; comp_image(:,:,1) = temp;
temp = comp_image(:,:,2); temp(cell_zimg(:,:,end)) = 0; comp_image(:,:,2) = temp;
temp = comp_image(:,:,3); temp(cell_zimg(:,:,end)) = 0; comp_image(:,:,3) = temp;


hit_mean = mean(cell_hit_prc);
hit_std = std(cell_hit_prc);
% profile off
% profile viewer
toc
%AVG_3 = mean(r);