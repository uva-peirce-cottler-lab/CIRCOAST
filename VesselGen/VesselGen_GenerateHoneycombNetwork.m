function [final_bw_vessel, final_bw_skel, stats] =VesselGen_GenerateHoneycombNetwork(img_dim, vess_rad_pix,umppix, ...
    restrict_rad,varargin)

% ARGUMENT PARSING
% What is not an Parameter is an optinal arg (param words must fail
% validation for add optional).
p = inputParser;
p.addParameter('SpawnPointNumLimit', [], @isnumeric);
p.addParameter('VesselLengthDensityLimits', [], @isnumeric);
p.addParameter('VesselAreaFractionLimits', [], @(x)(x>=0 && x<=1));
p.addParameter('MasterAngle', [], @(x)(x>=0 && x<=1));
p.addParameter('FigureVisible','on', @(x) strcmp(x,'off') || strcmp(x,'on'));
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

if isempty(SpawnPointNumLimit) && isempty(VesselAreaFractionLimits) &&...
        isempty(VesselLengthDensityLimits)
   error('Need to specify a limit (Vessel frac, length density, or spawn point number) for network');
end 

if VesselLengthDensityLimits(1)>50;
    [~,bw_skel_extra] = VesselGen_GenerateHoneycombNetwork(img_dim, vess_rad_pix,umppix, ...
    restrict_rad,'VesselLengthDensityLimits',[40 10]);
else
   bw_skel_extra=false(img_dim); 
end

isEven = @(x) mod(x,2)==0;
if isEven(restrict_rad); restrict_rad=restrict_rad+1; end

[bw_kernel, rk, ck] = VesselGen_RestrictKernel(restrict_rad,vess_rad_pix);

% Padded image coordinates
ri1 = 2*size(bw_kernel,1)+1;
ri2 = img_dim(1) +2*size(bw_kernel,1);
ci1 = 2*size(bw_kernel,2)+1;
ci2 = img_dim(2) +2*size(bw_kernel,2);

% Store points in 2D BW image
pad_pt_bw = false(img_dim + 2*size(bw_kernel)*2);
pad_ves_seg_bw = false(size(pad_pt_bw));
% Linear index of same points
pt_list=[];

% Possible spawn points in BW img 
pad_spawn_gs = zeros(size(pad_pt_bw),'uint8');
pad_spawn_gs(ri1-rk:ri2+rk, ci1-ck:ci2+ck)=1;
%Linear Indidicies of possible spawn points
static_pad_ind = 1:numel(pad_spawn_gs);
pad_spawn_ind = static_pad_ind(static_pad_ind);


% Loop until vesel density or fraction exceeded, or max points spawned
h = figure();
for n=1:1e3
    
    
    % Calculate new Spawn points
    pad_spawn_ind = static_pad_ind(logical(pad_spawn_gs));
    new_spawn_pt = pad_spawn_ind(randi([1 numel(pad_spawn_ind)], 1));
    
    [new_spawn_pt_rc(1) new_spawn_pt_rc(2)] = ...
        ind2sub(size(pad_spawn_gs),new_spawn_pt);
    
    % Add to list of spawn pts 
    pad_pt_ind(n) = new_spawn_pt; 
    pad_pt_rc(n,1:2) = new_spawn_pt_rc;
    
    % Add points, watershed [update pts, skel]
    pad_pt_bw(new_spawn_pt)=1;
    watershed_bw = ~watershed(bwdist(pad_pt_bw));
    bw_skel = bwmorph(bwmorph(watershed_bw(ri1:ri2,ci1:ci2),'Thin',Inf) | bw_skel_extra,'skel',Inf);
    
    % Remove pad with skel image
    npix_skel = sum(bw_skel(:));
    
    % Calculate vessel density
    calc_vd_mmpmm2(n) = npix_skel*(umppix/1000)/(img_dim(1)*umppix/1000).^2;
    
    
    % Estimate dilated area [update dilated vessels]
    est_vf(n) = (0.9*npix_skel*(2*vess_rad_pix+1) - n*pi*vess_rad_pix^2)/prod(img_dim);

    bw_vessel = bwdist(bw_skel) <=vess_rad_pix;
   
     
    calc_vf(n) = sum(bw_vessel(:))/numel(bw_vessel);
    
    % Restrict spawn indicies for next round
    pad_spawn_gs(pad_pt_rc(n,1)-(rk-1)/2:pad_pt_rc(n,1)+(rk-1)/2,...
        pad_pt_rc(n,2)-(ck-1)/2:pad_pt_rc(n,2)+(ck+1)/2-1) =  ...
        pad_spawn_gs(pad_pt_rc(n,1)-(rk-1)/2:pad_pt_rc(n,1)+(rk-1)/2,...
        pad_pt_rc(n,2)-(ck-1)/2:pad_pt_rc(n,2)+(ck+1)/2-1)- uint8(bw_kernel);
    % Generate new spawn indices
    
    if strcmp(FigureVisible,'on')
        figure(h);
        subplot(2,3,1); imshow(imadjust(pad_spawn_gs));
        
        subplot(2,3,2); imshow(imdilate(pad_pt_bw,strel('disk',3,0)));
        
        subplot(2,3,3); imshow(bw_skel);
        
        subplot(2,3,4); imshow(bw_vessel);
        
        subplot(2,3,5); plot(1:n,calc_vf(1:n));
        title({'Vessel Fraction',sprintf('%.3f',calc_vf(n))});
        hold on; plot(1:n,est_vf(1:n),'r'); hold off
        axis([0 n 0 0.5]);
        subplot(2,3,6); plot(1:n,calc_vd_mmpmm2(1:n));
        title({'Vessel Density',sprintf('%.2f',calc_vd_mmpmm2(n))});
        
        if strcmp(FigureVisible,'on'); pause(.01); end
    end
    
    % Test for loop exit conditions
    if  ~isempty(VesselLengthDensityLimits) && VesselLengthDensityLimits(1)< calc_vd_mmpmm2(n)
        fprintf('VesselLengthDensityLimits Reached.\n'); break;
    end
    
    if ~isempty(VesselAreaFractionLimits) && VesselAreaFractionLimits(1) < calc_vf(n)
        fprintf('VesselAreaFractionLimits Reached.\n'); break;
    end
    
    if ~isempty(SpawnPointNumLimit) && n >= SpawnPointNumLimit
       fprintf('SpawnPointNumLimit Reached.\n'); break 
    end
    
    % Ran out of spawn points
    if sum(logical(pad_spawn_gs))==0
        fprintf('No more spawn points available.\n'); break
        % Resize kernel
        % Recalculate spawn pts.
    end
 
end
close(h); 
% keyboard
% save('Vessel_regress_input.mat')
% load('Vessel_regress_input.mat')
% keyboard
if ~(isempty(VesselLengthDensityLimits) || VesselLengthDensityLimits(2)==0) ||...
        ~(isempty(VesselAreaFractionLimits) || VesselAreaFractionLimits(2)==0)
  [final_bw_vessel, final_bw_skel, rem_st] =...
    VesselGen_RegressNetwork(bw_skel, vess_rad_pix,umppix, ...
    'SpawnPointNumLimit',SpawnPointNumLimit,...
    'VesselLengthDensityTarget',diff(fliplr(VesselLengthDensityLimits)),...
    'VesselAreaFractionTarget',diff(fliplr(VesselAreaFractionLimits)));
else
    fprintf('No Initial regression on honeycomb network.\n')
    final_bw_vessel=bw_vessel;
    final_bw_skel = bw_skel;
    rem_st.vd_mmpmm2= calc_vd_mmpmm2(end);
    rem_st.vf =calc_vf(end); 
end

stats.VesselAreaFraction = rem_st.vf;
stats.VesselLengthDensity_mmpmm2 = rem_st.vd_mmpmm2;




% keyboard
end



% a=zeros(100,100);
% a(round(rand(1,30)*100^2))=1;
% d=bwdist(a);
% vn = imadjust(watershed(d));
% vn(logical(a))=1;
% imshow(vn)