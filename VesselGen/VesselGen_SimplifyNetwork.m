
function  bwf_skel = VesselGen_RemoveHiddenSegments(bw_skel);

img_dim = size(bw_skel);
    

%% Eliminate small holes in network
% Dilate all holes smaller than 20 pix and reskel
cc_holes = bwconncomp(~bw_skel,4);

bw_smhole = false(img_dim);
for n=1:cc_holes.NumObjects
    numel(cc_holes.PixelIdxList{n})
    if numel(cc_holes.PixelIdxList{n})<40
        bw_smhole(cc_holes.PixelIdxList{n})=1;
    end
end


%% Unite branchpoints close to one another
bw_bp = bwmorph(bw_skel,'branchpoints');
bp_strel = strel('disk',3,0);
cc_bp = bwconncomp(imdilate(bw_bp,bp_strel));
bw_close_bp = false(img_dim);
for n=1:cc_bp.NumObjects
%     numel(cc_bp.PixelIdxList{n})
    if numel(cc_bp.PixelIdxList{n})>sum(bp_strel.Neighborhood(:))
        bw_close_bp(cc_bp.PixelIdxList{n})=1;
    end
end



bw2_skel = bwmorph(bw_skel | imdilate(bw_smhole | bw_close_bp,strel('disk', 6,0)),'Thin',Inf);


bw3_skel = bwmorph(imdilate(bw2_skel,strel('disk',1,0)),'Thin',Inf) ;


%% Remove all endpoints
bw_border = false(img_dim); bw_border(1,:)=1; bw_border(:,1)=1;
bw_border(end,:)=1;bw_border(:,end)=1;
% bw_ep = bwmorph(bw3_skel,'endpoints') & bw_border;
bwf_skel=bw3_skel;
for n=1:100
    bw_rmv =  bwmorph(bwf_skel,'endpoints') & ~bw_border;
    if all(bw_rmv==0)
        break;
    end
    bwf_skel = bwf_skel & ~bw_rmv;
end


