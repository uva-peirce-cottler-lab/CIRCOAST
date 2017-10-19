function [bw_skel_out, complete_rem_seg_ind,rmv_seg_ind, MAX_SEGS_REMOVED] = ...
    VesselGen_TryRemoveSegment(bw_skel, rem_cc_seg, seg_rmv_ind,BW_Protected)


img_dim= size(bw_skel);
% bw_seg_temp=false(img_dim);
MAX_SEGS_REMOVED=0;
DBUG=0;


rem_bw_seg = false(img_dim);
for n=1:rem_cc_seg.NumObjects
    rem_bw_seg(rem_cc_seg.PixelIdxList{n})=1;
end

    if DBUG; figure; labelCCMatrix(rem_cc_seg,1:rem_cc_seg.NumObjects); end
    
% Scan through remove list, lowest score first, to find suitable
% segment to remove
for ind=1:rem_cc_seg.NumObjects
    
    % Temp BW for removing segment
    bw_skel_out = bw_skel;
    
    % Get indicies of segment to be removed
    seg_ind = seg_rmv_ind(ind);
    rem_pix_ind=rem_cc_seg.PixelIdxList{seg_ind};
    
     % Visualize segments proposed for removal
    complete_rem_seg_ind = VesselGen_CompleteNetworkSegment(bw_skel, rem_pix_ind);
    
    % Remove segment from skeleton
    bw_skel_out(complete_rem_seg_ind)=0;
%     bw_skel_out = bwmorph(bw_skel_out,'spur');
    
    bw_skel_out=bwareaopen(bw_skel_out,7);
    
    
    if DBUG; temp=false(img_dim); temp(rem_pix_ind)=1;
        imshow(im2uint8(cat(3,temp,rem_bw_seg, bw_skel_out))); end
    
    
    [VALID, pass_bv] = VesselGen_ValidNetwork(bw_skel_out);
        PROTECTED = any(BW_Protected(rem_pix_ind));

    
    % Print output
    if ind==1; fprintf('\tBPair,BTouch,CC_Rang\n');end
    fprintf('\t%i) SEG_%3.i:\t%i\t%i\t%i\t%i',ind,seg_ind,~PROTECTED,pass_bv(1),pass_bv(2),pass_bv(3));
    if PROTECTED; fprintf('_PROTECTED_'); end
    fprintf('\tAdding %i pixel to complete segment\n', ...
        numel(complete_rem_seg_ind)-numel(rem_pix_ind));


%     if  keyboard; end

    % Accept Changes
    if  VALID && ~PROTECTED
        fprintf('\tRemoving Seg: %i/%i, %i pixels\n',ind,...
            rem_cc_seg.NumObjects,numel(complete_rem_seg_ind));
        break
    elseif ind==rem_cc_seg.NumObjects
        fprintf('\tMax number of segments removed without violating rules.\n')
%         if rem_cc_seg.NumObjects>1; keyboard; end 
        % Release Protected Vessel
        [bw_skel_out, complete_rem_seg_ind,rmv_seg_ind, MAX_SEGS_REMOVED] = ...
            VesselGen_TryRemoveSegment(bw_skel, rem_cc_seg, seg_rmv_ind,false(img_dim));
        bw_skel_out=bwareaopen(bw_skel_out,7);
        if all(bw_skel_out==bw_skel); keyboard; end
        MAX_SEGS_REMOVED=1;
%         keyboard 
        break
    end
end
rmv_seg_ind = ind;

end