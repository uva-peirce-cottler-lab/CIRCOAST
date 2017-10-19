function [rmv_score] = VesselGen_ScoreNetworkSegments(rem_bw_skel, rem_cc_seg,score_arg_st)

% SMALL SCORE --> REMOVE SEGMENT
RemoveAngle=35;
    %TODO Protect longest couple of lines as VLD gets low. Algorithm?

fargs = fields(score_arg_st);
for n=1:numel(fargs); eval([fargs{n} '=' 'score_arg_st.' fargs{n} ';']);  end
img_dim = size(rem_bw_skel);

% Normalize score 0-100
scale100 = @(x) (x-min(x))*100/(max(x));
    
if SCORE_ANGLE
    orien_stats = regionprops(rem_cc_seg,'Orientation');
    seg_angle = arrayfun(@(x) x.Orientation,orien_stats);
    score_orien = ((seg_angle-RemoveAngle).^2)';
else
    score_orien=zeros(1,rem_cc_seg.NumObjects);
end




% Prevent long segments from dropping out
if SCORE_SEG_LEN
    seg_len = cellfun(@(x) numel(x), rem_cc_seg.PixelIdxList);
    score_seg_len = scale100(seg_len);
    
else
    score_seg_len = zeros(1,rem_cc_seg.NumObjects);
end


if SCORE_HOLE_SIZE
    % Score segment based on small the hole they border is first
    %   i.e. eliminate segments around small holes first
    cc_hole = bwconncomp(~rem_bw_skel,4);
    db_temp = zeros(img_dim);
    for n=1:cc_hole.NumObjects
        db_temp(cc_hole.PixelIdxList{n})=numel(cc_hole.PixelIdxList{n});
    end
    % Blend area values into vessel segments, want LOWER values to
    % supercede higher values, no have to reverse numbering for
    % dilation
    db_temp(db_temp==0)=max(db_temp(:))+1;
    rem_bw_hole_area = max(db_temp(:))-...
        imdilate(max(db_temp(:))-db_temp,strel('disk',2,0));
    seg_hole_area = cellfun(@(x) mean(rem_bw_hole_area(x)),rem_cc_seg.PixelIdxList);
    score_area = scale100(seg_hole_area);
    
else
    score_area=zeros(1,rem_cc_seg.NumObjects);
end


if SCORE_RAND
    score_rand = normrnd(50,10,[1 rem_cc_seg.NumObjects]);
else
    score_rand = zeros(1,rem_cc_seg.NumObjects);
end


if SCORE_ENDPOINT
    % Remove segments with endpoints first, 0 for adding EP (but not border points)
    bw_border=false(img_dim);
    bw_border(1,:) = 1; bw_border(:,1) = 1;
    bw_border(end,:)=1; bw_border(:,end) = 1;
    rem_bw_pad_ep = imdilate(bwmorph(rem_bw_skel,'endpoints')& ~bw_border,strel('disk',3));
    bw_seg_ep = rem_bw_pad_ep & rem_bw_skel;
    
    % If seg touches endpoint, give low score
    for n=1:rem_cc_seg.NumObjects
        seg_has_ep(n)= any(bw_seg_ep(rem_cc_seg.PixelIdxList{n})==1);
    end
    if max(seg_has_ep)==0
        score_endpoints =  zeros(1,rem_cc_seg.NumObjects);
    else
        score_endpoints=100-scale100(seg_has_ep);
    end
%     keyboard
else
   score_endpoints =  zeros(1,rem_cc_seg.NumObjects);
end

rmv_score = sqrt(score_orien.^2 + score_area.^2 +...
    4*score_endpoints.^2 + score_seg_len.^2 +score_rand.^2);
