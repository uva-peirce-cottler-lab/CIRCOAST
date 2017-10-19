function complete_seg_ind = VesselGen_CompleteNetworkSegment(bw_skel, rem_pix_ind)
img_dim=size(bw_skel);
img_ind = 1:prod(img_dim);
DBUG=0;

% Add missing penultimate endpoints to segment (had to pad the
% branchpoints to break up skeleton in 8 connected neighberhood)
bw_seg=false(img_dim);
bw_seg(rem_pix_ind)=1;
bw_dil_seg_temp = imdilate(bw_seg,strel('square',3));
% Find where padding overlaps with segment and endpoints to find
% pen ultimate enedpoint (but still must be on skeleton)
bw_bp = bwmorph(bw_skel,'branchpoints');
bw_dil_bp = imdilate(bw_bp,strel('square',3));
bw_ep = bwmorph(bw_skel,'endpoints');
bw_dil_ep = imdilate(bw_ep,strel('square',3));
bw_overlap = bw_dil_seg_temp & (bw_dil_bp | bw_dil_ep) & bw_skel;

if DBUG
temp=cat(3,bw_seg,bw_overlap,bw_ep);
figure;imshow(im2uint8(temp));
end


% Complete segment, find endpoint, trace pixels because CC does not trace
% pixel index (will need htis later is partially restoring this segment)
bw_complete_seg = bw_overlap | bw_seg;

% Include endpoint centers for segment, not pen ultimate
bw_complete_seg2 = (imdilate(bw_complete_seg,strel('square',3)) & bw_dil_ep & bw_skel) | bw_complete_seg; 

% Include endpoint centers for segment, not pen ultimate
bw_complete_seg3 = (imdilate(bw_complete_seg2,strel('square',3)) & bw_dil_ep & bw_skel) | bw_complete_seg2; 



if DBUG
temp=cat(3,bw_complete_seg3,bw_complete_seg3,bw_ep);
figure;imshow(im2uint8(temp));
end



bw_complete_seg_ep = bwmorph(bw_complete_seg3,'endpoints');
ep_ind = img_ind(bw_complete_seg_ep);

%See if overlap endpoints, extend to actual endpoints


%If no endpoints segment is a loop, choose overlap point next to
%branchpoint for starting point for trace
if isempty(ep_ind)
    ov_ind = img_ind(bw_overlap);
    ep_ind=ov_ind(1);
end


[r, c] = ind2sub(img_dim, ep_ind(1));

complete_seg_rc = bwtraceboundary(bw_complete_seg3,[r c],'N');
complete_seg_ind_2x = sub2ind(img_dim,complete_seg_rc(:,1),complete_seg_rc(:,2));

complete_seg_ind = complete_seg_ind_2x(1:sum(bw_complete_seg3(:)));

% Verifying its sill one CC
cc_complete_seg = bwconncomp(bw_complete_seg3);


% cc=bwconncomp(bw_skel);
% % temp=bw_skel;
% temp(cc_complete_seg);
% cc2=bwconncomp(temp);
% if cc
    


if cc_complete_seg.NumObjects~=1; keyboard; end
if DBUG
temp=cat(3, bw_seg, bw_complete_seg3,bw_bp);
figure;imshow(im2uint8(temp));
fprintf('Adding %i pixels, %i objects to segment.\n',...
    numel(cc_complete_seg.PixelIdxList{1}),cc_complete_seg.NumObjects);
end
% keyboard
