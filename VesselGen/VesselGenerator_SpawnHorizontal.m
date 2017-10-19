function bw_vessels = VesselGenerator_SpawnHorizontal(varargin)


% EXAMPLE: 
%   bw_vessels = VesselGenerator_SpawnHorizontal([1024 1024], 20, 2, 'VesselDensity', .1)


%% ARGUMENT PARSING
p = inputParser;
p.addRequired('img_size', @(x) numel(x) == 2)
p.addRequired('vessel_diam_um', @(x) numel(x) == 1 && isnumeric(x))
p.addRequired('umppix', @(x) numel(x) == 1 && isnumeric(x))
p.addParamValue('NumVessels', 1000, @isnumeric);
p.addParamValue('VesselDensityFraction', 10, @isnumeric);
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n = 1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% Calculate vessel diameter in pixels
vessel_diam_pix =  round(vessel_diam_um/umppix);

% Index to record porbability of spawning new blood vessel
spawn_ind_prob = ones(1, img_size(1));

% Temp image to generate current vessel, add to bw_vessels each iteration
temp = false(img_size);
bw_vessels = false(img_size);
bw_skels = false(img_size);

% Record vector of vessel length and cumulative network length
vessel_length = zeros(1, NumVessels);
vessels_area = 0;

% Kernel used to move
move_kernel = [0 0 1; 0 0 1; 0 0 1];

for n = 1 : NumVessels
    % Generate renadom scores, multiply by spawn_ind, previous vessel spawn
    % points reduce the surrounding spawn+ind
    [~, ri] = max(rand([1 numel(spawn_ind_prob)]).*spawn_ind_prob);
    
    % Initiate Vectors for traversing the image
    rv = ri; cv = 1; 
    for k = 1:(512*5)
        bounds_kernel = [rv(k)>1 rv(k)>1 rv(k)>1; ...
            0 0 0;
            rv(k)<img_size(1) rv(k)<img_size(1) rv(k)<img_size(1)];
        [rshift, cshift] = VesselGenerator_CalculateMove(move_kernel .* bounds_kernel);
        rv(k+1) = rv(k) + rshift;
        cv(k+1) = cv(k) + cshift;
        
        if rv(end) == 0 || cv(end) == 0 || ...
                rv(end) == (img_size(1) + 1) || cv(end) == (img_size(2) + 1);
            break;
        end
    end
    rv(end) = []; cv(end) = [];
    temp(sub2ind(img_size, rv,cv)) = 1;
    temp_dil = imdilate(temp,strel('disk',ceil((vessel_diam_pix+1)/2),0));
    
    % Record locatin where vessel spawned (2x its thickness)
    spawn_rows = imdilate(temp_dil(:,1), ones(vessel_diam_pix,1));
    
    % Subsequent vessels are 1/2as likely to spawn at the current location
    spawn_ind_prob(spawn_rows) = spawn_ind_prob(spawn_rows)/3;

    % Add vessel to image of all vessels, reset temp image
    bw_skels = bw_skels | temp;
    bw_vessels = bw_vessels | temp_dil;
    temp(:)=0;
    
    % Record vessel area
    vessels_area = sum(bw_vessels(:));
    
    % Exit if threshold vessel density reached
    if vessels_area/prod(img_size) >= VesselDensityFraction
        fprintf('Vessel Density reached after %0.f vessels\n', n);
        fprintf('Target Vessel Density: %.4f\n', VesselDensityFraction);
        fprintf('Actual Vessel Density: %.4f\n', vessels_area/(prod(img_size)));
       break 
    end
  
end


