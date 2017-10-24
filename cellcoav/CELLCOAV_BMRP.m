function out_st = CELLCOAV_BMRP(varargin)
% CELLCOAV_BMRP(imgs_cell, cell_diam_um, umppix, n_obs_success,tot_obs_trials)
%% ARGUMENT PARSING
p = inputParser;
p.addRequired('imgs_cell', @(x) iscell(x) || isnumeric(x) || ...
    islogical(x) || isempty(x));
p.addRequired('cell_diam_um', @isnumeric);
p.addRequired('umppix', @isnumeric);
p.addRequired('tot_obs_trials', @isnumeric);
p.addParamValue('n_obs_success',[], @isnumeric);
p.addParamValue('p_success', [],@isnumeric);
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% Check to make sure that sufficient data is specified
assert(~isempty(imgs_cell) || ~isempty(p_success),...
    'Either imgs_cell or p_success must be specified (nonzero)');

% If imgs_cell is an img, make it a cell since we expect a cell of images
if (isnumeric(imgs_cell) || islogical(imgs_cell)) && all(size(imgs_cell>1)); 
    imgs_cell = {imgs_cell}; 
end


% Calculate success rate across image, assume images may be diff. pixel dimensions
% if isempty(p_success)
p_success = CELLCOAV_CellDilateVessFrac(imgs_cell,cell_diam_um,umppix);
% fprintf('BMRP ICF:%.3f\n',p_success)
% end

% General Method Calculating Mean from PDF
bino_pdf = binopdf(0:tot_obs_trials, ...
    tot_obs_trials,p_success);
bino_cdf = binocdf(0:tot_obs_trials, ...
    tot_obs_trials,p_success);

if ~isempty(n_obs_success)
    % Two tail binomial hypothesis test: greater than random
    % p = 2*(1-cdf(bino)),   +1 for 0 based indexing in cdf    

    % One tailed test, greater than random
    out_st.bmrp_cellcoav_p = min([(1-bino_cdf(round(n_obs_success+1)))]);% ...
end

% Long way of calculating binomial mean
% bino_mean = sum(discr_bino_prob .* (0:ntrials))/...
%     ntrials;

% Minimum number of trials for statiscial significance 
% to determine whether result is less than or greater than random
for n = 1:500
    discr_bino_cdf = binocdf(0:n, n,p_success);
    if any(discr_bino_cdf>.975 & discr_bino_cdf<1) && ...
        any(discr_bino_cdf<.025 & discr_bino_cdf>0); break; end
end
binom_min_ntrials = n;

% Binomail Distribution mean and std
%         u = n*p/n
binom_frac_mean = p_success;
%         std = sqrt(n*p*(1-p))
binom_frac_std = sqrt(tot_obs_trials*p_success*(1-p_success))/tot_obs_trials;


out_st.bmrp_icf_mean = binom_frac_mean;
out_st.bmrp_icf_std = binom_frac_std;
out_st.bmrp_min_ntrials = binom_min_ntrials;



