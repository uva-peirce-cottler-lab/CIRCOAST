function [VALID, pass_bv] = VesselGen_ValidNetwork(bw_skel_out)

img_dim =size(bw_skel_out);
 
% Find new edge points
bw_border=false(img_dim); 
bw_border(1,:)=1;bw_border(end,:)=1;
bw_border(:,1)=1; bw_border(:,end)=1;
bw_edgepoints = bwmorph(bw_skel_out,'endpoints') & bw_border;

% Find segments
bw_segs = bw_skel_out&~imdilate(bwmorph(bw_skel_out,'Branchpoints'),strel('square',3));
cc_seg = bwconncomp(bw_segs);

% Check that cc has not grown to much
cc_skel_temp = bwconncomp(bw_skel_out);

% Each CC must pass through border of image
n_solo_ep_pairs = zeros(1,cc_skel_temp.NumObjects);
n_edge_ep_pts = zeros(1,cc_skel_temp.NumObjects);
for m = 1:cc_skel_temp.NumObjects
    %             must have two endpoints greater than 1/2 distance away from each other
    cc_ind = cc_skel_temp.PixelIdxList{m};
    bv = bw_edgepoints(cc_ind);
    [r, c] = ind2sub(img_dim,cc_ind(bv));
    %             temp=false(img_dim); temp(cc_ind(bv))=1; imshow(temp)
    pix_dist = pdist([r c]);
    n_edge_ep_pts(m) = size([r c],1);
    n_solo_ep_pairs(m) = sum(pix_dist > img_dim(1)/1.5);
end

% Atleast 1 or 2 CCs must touch endpoints
CC_CONNECT_BORDER_PASS = (sum(n_solo_ep_pairs>=1) > 0) || (cc_seg.NumObjects<=1 && cc_skel_temp.NumObjects<=1);
% cc_seg.NumObjects
% cc_skel_temp.NumObjects

CC_TOUCH_BORDER_PASS = sum(n_edge_ep_pts>=1)==cc_skel_temp.NumObjects;
ENDPOINT_NUM_PASS =1;%Atleast 2 CCs
CC_RANGE_PASS = cc_skel_temp.NumObjects>=0 && cc_skel_temp.NumObjects<=3;
%         keyboard
%

pass_bv = [CC_CONNECT_BORDER_PASS CC_TOUCH_BORDER_PASS CC_RANGE_PASS];

VALID = all(pass_bv);


% keyboard
end