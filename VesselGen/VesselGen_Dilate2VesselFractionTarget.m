function [bw_vessel_frac,obs_vessel_frac] = ...
    VesselGen_Dilate2VesselFractionTarget(bw, vessel_frac,cell_diam_um,...
    umppix,CELL_AREA_FRACT)
DBUG=0;

if CELL_AREA_FRACT
    cellAreaFract = @(x) CELLCOAV_CellDilateVessFrac(x,cell_diam_um,umppix);
else
    cellAreaFract = @(x) sum(x(:))/numel(x);
end

% Euclidean distance transform
ed_gs = bwdist(bw);

% Initially dilate image until it exceed vessel frac, then step back once
for n=0:512
   %Threshold euclidean distance from vessel
   bw = ed_gs <= n;
   vf(n+1) = cellAreaFract(bw);
   if vf(n+1) > vessel_frac; break; end 
end
assert(n~=0,sprintf('Vessel Image already exceed vessel frac target: %.2f<%.2f',vf(n+1),vessel_frac))
% Initial Vessel Image
bw_vessel_frac_rough = ed_gs <= n-1;

% Calculate how many pix needed
tot_pix_required = round(numel(bw).*vessel_frac);

% Identify boundaries of vessels to partially fill
bw_border = bwdist(bw_vessel_frac_rough)==1;
B_trace = bwboundaries(~bw_vessel_frac_rough);
% keyboard


% Clear out tracings on border of image
for n=1:numel(B_trace)
    ind= sub2ind(size(bw),B_trace{n}(:,1),B_trace{n}(:,2));
    B_border{n}=B_trace{n}(bw_border(ind),:);
end


bw_vessel_frac = bw_vessel_frac_rough;
for n=1:numel(B_border)
    prev_npix(n) = cellAreaFract(bw_vessel_frac).*numel(bw);
   
    % Convert RC coordinates to linear indicies of border
    border_ind = sub2ind(size(bw),B_border{n}(:,1),B_border{n}(:,2));
    
    % Try adding the border to a test image
    bw_test = bw_vessel_frac;
    bw_test(border_ind)=1;
    
    %Calculate pixels after dilation
    curr_npix(n) = cellAreaFract(bw_test).*numel(bw);
    diff_npix(n) = curr_npix(n) - prev_npix(n);
    diff_target_npix(n) = tot_pix_required-curr_npix(n);
    
    
    
    % If less pixels needed than filling in border CC, pertially fill
    if diff_target_npix(n) <=0
        fprintf('Partially restoring border %i.\n',n);
        % Loop here to find exact number of pixels need to be added so
        % after dilation meets vessel frac target
        extra_npix = ceil(numel(border_ind).*...
            (1-(-diff_target_npix(n) / (curr_npix(n) -prev_npix(n)))) );
        for k=1:1e3
            % Fill in extra border pixels
            bv = false(size(border_ind));
            bv(1:extra_npix(k))=1;
            % Create new bw with added pixels
            bw_test = bw_vessel_frac;
            bw_test(border_ind(bv))=1;
            
            % Calc remaining pixels
            npix_remain(k) = tot_pix_required - cellAreaFract(bw_test)*numel(bw_test);
            % If no pixels remaining, then exit loop
            if npix_remain(k)==0; break; end
            
            % Iterate extra pixels needed convervatively
            extra_npix(k+1) = extra_npix(k) + sign(npix_remain(k));
            if k>2 && extra_npix(k+1)==extra_npix(k-1); 
                fprintf('Oscillating values\n')
                break; 
            end
           
        end 
        % Commit Changes
        bw_vessel_frac= bw_test;
        fprintf('\t%i Pixels added to border region %i to reach VF\n',extra_npix(k),k)
        % BW is now adjsuted for area fraction, exit for loop 
        break;
    else
      % Commit changes
      bw_vessel_frac=bw_test;
    end
    
end
% Calculate actual VF to verify
obs_vessel_frac = cellAreaFract(bw_vessel_frac);

if DBUG
    imshow(im2uint8(cat(3,bw_vessel_frac&~bw,bw,false(size(bw)))))
end

% Test if vessel fraction is to closest pixel
assert(abs(obs_vessel_frac-vessel_frac)<1e3,'Failed to reach vessel fraction target');

% keyboard

end