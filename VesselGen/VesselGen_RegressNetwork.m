function [rem_bw_vessel, rem_bw_skel, stat_st] =...
    VesselGen_RegressNetwork(bw_skel, vess_rad_pix, umppix, varargin)
DBUG=0;
img_dim =size(bw_skel);

if isempty(bw_skel); error('VesselGen_RegressNetwork: input bw_skel is empty'); end
% ARGUMENT PARSING
% What is not an Parameter is an optinal arg (param words must fail
% validation for add optional).
p = inputParser;
p.addParameter('SpawnPointNumLimit', [], @(x) isnumeric(x) || isempty(x));
p.addParameter('VesselLengthDensityTarget', [], @(x) isnumeric(x) || isempty(x));
p.addParameter('VesselAreaFractionTarget', [], @(x) all(x>=0 & x<=1) || isempty(x));
p.addParameter('Cell_Rad_Pix', 0,@isnumeric)
p.addParameter('RemoveAngle', [], @(x)(x>=0 && x<=1));
p.addParameter('EXACT_SEG_RESTORE', 1, @(x)(x==0 || x==1));
p.addParameter('BW_Protected', false(size(bw_skel)), @(x) islogical(x) || numel(x)==numel(bw_skel));
p.addParameter('SCORE_HOLE_SIZE', 1, @(x)(x==0 || x==1));
p.addParameter('SCORE_ANGLE', 1, @(x)(x==0 || x==1));
p.addParameter('SCORE_RAND', 1, @(x)(x==0 || x==1));
p.addParameter('SCORE_SEG_LEN', 0, @(x)(x==0 || x==1));
p.addParameter('SCORE_ENDPOINT', 1, @(x)(x==0 || x==1));
p.addParameter('FigureVisible','on', @(x) strcmp(x,'off') || strcmp(x,'on'));
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

for n=find(cellfun(@(x) ~isempty(regexp(x,'SCORE_','once')),fargs))'
    score_arg_st.(fargs{n})=p.Results.(fargs{n});
end
MAX_SEGS_REMOVED=0;

% Define function to calculate vessel frac (with or without cell dilation.
imDilFcn=@(x) bwdist(bwdist(x) <= vess_rad_pix) <= Cell_Rad_Pix;


% Check that # CC's alright for skel
cc_check = bwconncomp(bw_skel);
[~, ix]=sort(cellfun(@(x) numel(x),cc_check.PixelIdxList),'ascend');
for n=1:cc_check.NumObjects-3
    bw_skel(cc_check.PixelIdxList{ix(n)})=0;
end

% Keep record of original segments and endpoints
bw_orig_bp = bwmorph(bw_skel,'branchpoints');
% Edgepoints are endpoints onthe border of the image
bw_border=false(img_dim); bw_border(1,:)=1;bw_border(end,:)=1;
bw_border(:,1)=1; bw_border(end,:)=1;
bw_orig_edgepoints = bwmorph(bw_skel,'endpoints')& bw_border;
bw_orig_seg= bw_skel & ~imdilate(bw_orig_bp,strel('square',3));


if DBUG; imshow(im2uint8(cat(3, bw_orig_bp,bw_orig_seg,bw_skel))); end
% cc_orig_seg = bwconncomp(bw_orig_seg);

% Initialize skel and remember original endpoints (these go to img border)
rem_bw_skel = bw_skel;
rem_cc_skel = bwconncomp(rem_bw_skel);
bw_seg_temp=false(img_dim);

hf = figure('Visible',FigureVisible);
for k=1:1e6
    fprintf('Round %i:\n',k);
    rem_bw_skel = bwmorph(rem_bw_skel,'skel',Inf);
    if isempty(rem_bw_skel); keyboard; end
    
    % Recalc brnachpoints (they change as segments removed
    rem_bw_bp = bwmorph(rem_bw_skel,'branchpoints');
    rem_bw_seg = rem_bw_skel & ~imdilate(rem_bw_bp,strel('square',3));
    cc_seg_pre_rmv = bwconncomp(rem_bw_seg);
    
    if DBUG; imshow(im2uint8(cat(3,rem_bw_skel,rem_bw_seg, rem_bw_bp))); end
    
     
    % Find Endpoints
    bw_border=false(img_dim); bw_border(:,1)=1; bw_border(1,:)=1;
    bw_border(:,end)=1; bw_border(end,:)=1;
    bw_endpoints = bwmorph(rem_bw_skel,'endpoints') & ~bw_border;
    bw_edgepoints = bwmorph(rem_bw_skel,'endpoints') & bw_border;
    
    % Score segments on whether they should be removed
    [rmv_score] = VesselGen_ScoreNetworkSegments(rem_bw_skel, cc_seg_pre_rmv,score_arg_st);
    [srt_rmv_score, srt_rmv_ix] = sort(rmv_score);
    
    % Visualize segments proposed for removal
    if DBUG; figure; labelCCMatrix(cc_seg_pre_rmv,1:cc_seg_pre_rmv.NumObjects); end
    
    % Check validity of network 
    [VALID, pass_bv] = VesselGen_ValidNetwork(rem_bw_skel);
    if ~VALID; keyboard;    
    temp=false(img_dim); temp(prev_rem_seg_ind)=1;
    imshow(im2uint8(cat(3,temp,rem_bw_skel,rem_bw_skel)));
    end
    
    bw_skel_pre_rmv = rem_bw_skel;
    cc_skel_pre_rmv = bwconncomp(rem_bw_skel);
    
    [rem_bw_skel,rem_seg_ind,rmv_seg_ind, MAX_SEGS_REMOVED] = VesselGen_TryRemoveSegment(...
        bw_skel_pre_rmv, cc_seg_pre_rmv, srt_rmv_ix, BW_Protected);

    if DBUG   
        temp=false(img_dim); temp(cc_seg_pre_rmv.PixelIdxList{rmv_seg_ind})=1;
        rgb_temp=cat(3,temp,rem_bw_skel, temp);
        figure;imshow(im2uint8(rgb_temp));
    end
    % Get CC # After removal
    cc_post_rmv = bwconncomp(rem_bw_skel);
    
    
    % UPDATE vessel image from skeleton
    rem_bw_vessel = bwdist(bwdist(rem_bw_skel) <= vess_rad_pix);
    
    %Calculate metrics
    calc_vf(k) = sum(sum(imDilFcn(rem_bw_skel)))./prod(img_dim);
    
    % VLD: pix * mm_per_pix/area_of_image_mm2
    %                        npix_skel*(umppix/1000)/(img_dim(1)*umppix/1000).^2;
    calc_vd_mmpmm2(k) =  sum(rem_bw_skel(:))*(umppix/1000)/(img_dim(1)*umppix/1000).^2;
    fprintf('\tVLD: %4.3f, VF:%0.2f\n',calc_vd_mmpmm2(k),calc_vf(k));
    
    
    %Display Results
    if strcmp(FigureVisible,'on')
        temp=false(img_dim); temp(rem_seg_ind)=1;
        rgb_visualize = im2uint8(cat(3,temp,rem_bw_skel,...
            BW_Protected));
        figure(hf); imshow(rgb_visualize)
        %     subplot(2,2,1); imshow(rem_bw_skel);
        %     subplot(2,2,2); %imshow(rgb_visualize);
        %     subplot(2,2,3); plot(1:k,rem_calc_vf(1:k));
        %     title({'Vessel Fraction',sprintf('%.3f',rem_calc_vf(k))});
        %     axis([0 n 0 0.5]);
        %     subplot(2,2,4); plot(1:k,rem_calc_vd_mmpmm2(1:k));
        %     title({'Vessel Density',sprintf('%.2f',rem_calc_vd_mmpmm2(k))});
        pause(.1);
    end
    %
    
    
    % For debugging pruporses, rememebr previous segment removed
    prev_rem_seg_ind = rem_seg_ind;

    %% Test for loop exit conditions
    if  ~isempty(VesselLengthDensityTarget) && ...
            (VesselLengthDensityTarget> calc_vd_mmpmm2(k) || MAX_SEGS_REMOVED)
        fprintf('VesselLengthDensityTarget Reached.\n');
%         keyboard
        % Recreate vesssel image
        rem_bw_vessel= bwdist(bwdist(rem_bw_skel) <= vess_rad_pix);
        if EXACT_SEG_RESTORE
            fprintf('Partially restoring last removed segment for precision.\n');
            % Fill last segment in partially to reach goal to the pixel, recalc imgs and metrics
            % mm/mm2 -> vld*1000/ummpix /(ummpix*img_dim(1)/1000).^2
            vld_mmpmm2_to_pix  = @(vld) vld*1000/umppix * (umppix*img_dim(1)/1000).^2;
            
            % Calculate remaining pixels needed
            rem_pix = round(vld_mmpmm2_to_pix(VesselLengthDensityTarget)-...
                vld_mmpmm2_to_pix(calc_vd_mmpmm2(k)));
            
            % Try restoring one side of removed segment
            bw_restore_skel = rem_bw_skel;
            restore_seg_ind = rem_seg_ind;
            % Sometimes pixels removed in bwareaopen
            bw_restore_skel(restore_seg_ind(1:min([rem_pix numel(restore_seg_ind)])))=1;
            
            
            % For correct restorating of segment, network must be valid and
            % same number of CC's compared to b/f segment removal
            cc_post_restore = bwconncomp(bw_restore_skel);
            [VALID, pass_bv] = VesselGen_ValidNetwork(bw_restore_skel);

%             imshow(im2uint8(cat(3, bw_restore_skel, rem_bw_skel,false(img_dim))))
            if (cc_post_restore.NumObjects ~= cc_skel_pre_rmv.NumObjects) ...
                    || ~VALID
                bw_restore_skel = rem_bw_skel;
                restore_seg_ind = flipud(rem_seg_ind);
                bw_restore_skel(restore_seg_ind(1:rem_pix))=1;
                cc_post_restore = bwconncomp(bw_restore_skel);
                [VALID, pass_bv] = VesselGen_ValidNetwork(bw_restore_skel);
                
                if ~VALID; 
%                     temp=false(img_dim); temp(restore_seg_ind(1:rem_pix))=1;
%                     figure;
%                     imshow(im2uint8(cat(3, rem_bw_skel, bw_restore_skel, temp)));
%                     keyboard; 
                end
            end
            
            
            
            if cc_post_restore.NumObjects ~= cc_skel_pre_rmv.NumObjects
%                 keyboard; 
                temp=false(img_dim); temp(rem_seg_ind(1:rem_pix))=1;
                temp2 = false(img_dim); temp2(rem_seg_ind)=1;
                rgb_temp=cat(3,temp, rem_bw_skel,temp2);
%                 figure;imshow(im2uint8(rgb_temp));
            end 
            
            % Commit changes
            rem_bw_skel = bw_restore_skel;
            % Recreate vesssel image
            rem_bw_vessel=bwdist(rem_bw_skel) <= vess_rad_pix; 
            
        end
        break;
    end
    
    if ~isempty(VesselAreaFractionTarget) && ...
            (VesselAreaFractionTarget > calc_vf(k) || MAX_SEGS_REMOVED)
        fprintf('VesselAreaFractionTarget Reached: %.3f.\n',VesselAreaFractionTarget);
        if EXACT_SEG_RESTORE
            bw_dil = imDilFcn(rem_bw_vessel);
            npix_diff = (VesselAreaFractionTarget-calc_vf(k))*prod(img_dim);
            dil_rad = Cell_Rad_Pix+vess_rad_pix;
            
           
            % Restoration indicies are last removed
            restore_seg_ind = rem_seg_ind;
            bw_restore_skel = rem_bw_skel;
            % Try restoring half os segment see if network still valid
            bw_restore_skel(restore_seg_ind(ceil(numel(restore_seg_ind)/2)))=1;
            % If not valid, then flip restore indices
            [VALID, ~] = VesselGen_ValidNetwork(bw_restore_skel);
            if ~VALID; restore_seg_ind=restore_seg_ind(fliplr(1:numel(restore_seg_ind)));
            end
            
            %if not split restore segment indices
            skel_npix_diff = floor(npix_diff./(2*dil_rad+1));
            restore_vf=[];
            for kk = 1:numel(restore_seg_ind)*2
%                 fprintf('Restore: %i pixels.\n',skel_npix_diff(kk))
                if kk>3 && (skel_npix_diff(kk)==skel_npix_diff(kk-1) ||...
                        skel_npix_diff(kk)==skel_npix_diff(kk-2))
                    fprintf('\tRestoring Segment oscillating, exiting...\n')
                    break
                end
                
                % Try restoring one side of removed segment
                bw_restore_skel = rem_bw_skel;
                
                % Restore part of segment
                bw_restore_skel(restore_seg_ind(...
                    1:min([skel_npix_diff(kk) numel(restore_seg_ind)])))=1;
%                  bw_restore_skel(restore_seg_ind)=1;
                 
                bw_restore_vessel = imDilFcn(bw_restore_skel);
                restore_vf(kk) = sum(bw_restore_vessel(:))./numel(bw_restore_vessel);
           
                skel_npix_diff(kk+1)= skel_npix_diff(kk)+...
                    sign(VesselAreaFractionTarget-restore_vf(kk));
                
            end
            
            [VALID, pass_bv] = VesselGen_ValidNetwork(bw_restore_skel);
%             if ~VALID; keyboard;end
           
            keyboard
            fprintf('Partially restoring last removed segment for precision.\n');
            % Fill last segment in partially to reach goal to the pixel
            % Fill last segment in partially to reach goal to the pixel, recalc imgs and metrics
            rem_bw_skel =bw_restore_skel;
            rem_bw_vessel= bwdist(bw_restore_skel) <= vess_rad_pix;
        end
        break;
    end
    %     keyboard
end
close(hf);

% Recalc metrics in case of partial refilled segment (EXACT_SEG_RESTORE)
calc_vf(k+1) = sum(rem_bw_vessel(:))/numel(rem_bw_vessel);
calc_vd_mmpmm2(k+1) =  sum(rem_bw_skel(:))*(umppix/1000)/(img_dim(1)*umppix/1000).^2;


stat_st.vf=calc_vf(end);
stat_st.vld_mmpmm2=calc_vd_mmpmm2(end);



end