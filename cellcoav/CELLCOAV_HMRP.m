function out_st = CELLCOAV_HMRP(varargin)
p = inputParser;
p.addRequired('imgs_cell', @(x) iscell(x) || isnumeric(x) || ...
    islogical(x) || isempty(x));
p.addRequired('cell_diam_um', @isnumeric);
p.addRequired('umppix', @isnumeric);
p.addRequired('tot_obs_trials', @isnumeric);
p.addParameter('n_obs_success',[], @isnumeric);
p.parse(varargin{:});
% Import parsed variables into 
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

if ~iscell(imgs_cell);
imgs_cell = {imgs_cell};    
end
% Closest odd number
closestOdd = @(x) round(round((x-1)/2))*2+1;

cell_diam_pix = closestOdd(cell_diam_um/umppix);


% Load hypergeomtric table for packing ratios
hgy_pack_rat_tbl = getappdata(0,'hgy_pack_rat_tbl');
if isempty(hgy_pack_rat_tbl)
    hgy_pack_rat_tbl = readtable([getappdata(0,'proj_path') '/ArcasGui/hyg_pack_rat.csv']);
    setappdata(0,'hgy_pack_rat_tbl',hgy_pack_rat_tbl);
end
% Look up packing ratio value
ind = find(cell_diam_pix==hgy_pack_rat_tbl.cell_diam_pix);
if cell_diam_pix> max(hgy_pack_rat_tbl.cell_diam_pix);
    ind = size(hgy_pack_rat_tbl.cell_diam_pix,1);
elseif isempty(ind)
    error('HMRP Kernel size not in table');
end
 
pack_rat = hgy_pack_rat_tbl.mean_pack_rat(ind);
assert(cell_diam_pix==hgy_pack_rat_tbl.cell_diam_pix(ind) ||...
    cell_diam_pix>hgy_pack_rat_tbl.cell_diam_pix(end),...
    'Error: cell diam (pix) did not look up from tbl correctly');

% Calculate success rate across image, assume images may be diff. pixel
% dimensions
[p_success, npix_tot] = CELLCOAV_CellDilateVessFrac(imgs_cell,cell_diam_um,umppix);


% hypgcdf(x,N,k,n)
% computes the hypergeometric cdf at each of the values in 
%  N: Population Size, max theoretical # of cells in image
N =  round(sum(npix_tot)*pack_rat/(pi*cell_diam_pix));  
%  k: # cells with the desired characteristic in the populatio
k = round(N*p_success);
%  n: and number of samples drawn
n = tot_obs_trials;
%  x: # of items in sample classified as successes
hypg_pdf = hygepdf(0:n,N,k,n);
hypg_cdf = cumsum(hypg_pdf);


% h(x; N, n, k) = [ kCx ] [ N-kCn-x ] / [ NCn ]
% http://stattrek.com/probability-distributions/hypergeometric.aspx
%  N: The number of items in the population.
%  k: The number of items in the population that are classified as successes.
%  n: The number of items in the sample.
%  x: The number of items in the sample that are classified as successes.
% The mean of the distribution is equal to n * k / N .
hmrp_icf_mean = n * k / (N*n);

% STD of Hypergeometric
% The variance is n * k * ( N - k ) * ( N - n ) / [ N2 * ( N - 1 ) ] 
hmrp_icf_std = sqrt(n * k * ( N - k ) * ( N - n ) / (N.^2 * ( N - 1 ) ))/n;

% Two tail binomial hypothesis test: greater than random
% p = 2*(1-cdf(bino)),   +1 for 1 based indexing
if ~isempty(n_obs_success)
    out_st.hmrp_cellcoav_p = (1-hypg_cdf(n_obs_success+1));
end

% Minimum number of trials to determine whether result is less than or
% greater than random
for m = 1:500
    discr_hypg_cdf = hygecdf(0:tot_obs_trials,k,m,N);
    if any(discr_hypg_cdf>.95 & discr_hypg_cdf<1) && ...
        any(discr_hypg_cdf<.05 & bino_cdf>0); break; end
end
hypg_min_ntrials = m; 


out_st.hmrp_icf_mean = hmrp_icf_mean;
out_st.hmrp_icf_std = hmrp_icf_std;
% out_st.hypg_min_ntrials = hypg_min_ntrials;

fprintf('Hyg mean ICF:%.4f\n',hmrp_icf_mean);
% keyboard
