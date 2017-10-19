function [live_mean, live_std, dead_mean, dead_std] = ...
ArcasGui_RandomExperimentSimulator(img_dir, csv_name, csv_path)
%UNTITLED Summary of this function goes here
%   This simulation acts as a true negative for a cell association
%   experiment. It is basically  simulated control, evaluating what the
%   colocalization score would be if the experimental data was just one of
%   the monte carlo trials

if ~exist('img_name','var');
    [csv_name, csv_path] = uigetfile('*.csv','Select CSV of Manually Counted Results');
    img_dir = uigetdir(pwd,'Select Folder of Images');
end

%% Parse CSV file containing experimental data
csv_txt = fileread([csv_path csv_name]);
csv_lines = regexp(csv_txt, '\r\n|\n\r|\n','split');
csv_fields = cellfun(@(x) regexp(x,',','split'),csv_lines,'UniformOutput',0);
nfieldnames = 4;
csv_fields(cellfun(@(x) numel(x) ~= nfieldnames,csv_fields))=[];

csv_data = vertcat(csv_fields{:});

for c = 1:size(csv_data,2);
    if isnan(str2double(csv_data{2,c}))
        exp_st.(csv_data{1,c}) = csv_data(2:end,c);
    else exp_st.(csv_data{1,c}) = cellfun(@(x) str2double(x),csv_data(2:end,c));
    end
    
end

n_imgs = numel(exp_st.img_name);


umppix = 1.2656;
% event_diam_um = 8;
% event_diam_um = 11;
tot_trials = 10000;
nboots = tot_trials;%1000;
n_resample = 10000;
max_trials_per_cycle = 1500;
handles = [];



% keyboard
for n = 1:n_imgs
    
    img = imread([img_dir '/' exp_st.img_name{n} '.tif']);
    blue_bw = img(:,:,3) > 0;
    
    % Change cell diameter based on filename
    temp = regexp(exp_st.img_name{n},'.*?_.*?(.)_','tokens','once');
    isLive = temp{1} == 'L';
    if isLive; cell_diam_um = 11; else cell_diam_um = 8; end
    events_per_trial = exp_st.total(n);
    
    fprintf('%s: isLive: %.f, events: %0.f\n', exp_st.img_name{n}, isLive,events_per_trial );

    % Record Img name, whether live/dead
    sim_st.img_name{n} = exp_st.img_name{n};
    sim_st.isLive(n) = isLive;
    sim_st.exp_hits(n) = exp_st.hits(n);
    sim_st.exp_misses(n) = exp_st.misses(n);
    sim_st.exp_total(n) = exp_st.total(n);
    
    
    % Calculate mean and std with binomial distribution
     [sim_st.binom_frac_means(n), sim_st.binom_frac_stds(n)] = ...
           CELLCOAV_BMRP(blue_bw, cell_diam_um, umppix,exp_st.total(n));
    
    
    % Randomly generate 1000 binomial random fractions, 
    brf = binornd(exp_st.total(n),sim_st.binom_frac_means(n), [1 tot_trials])/exp_st.total(n);
    ind = randi(tot_trials,nboots);
    bino_mean = mean(brf);
    bino_std = std(brf);
    
    % Colocalizato score: (obs coloc frac) / mean(binom coloc frac)
    % In this negatice control case:
    %       score = (single binom coloc fraction) / mean(binom coloc frac)
    sim_st.coloc_score_rand(1:nboots,n) = (brf ./ mean(brf))';
    
end
n_dead = sum(~sim_st.isLive);
n_live = sum(sim_st.isLive);

% Seperate dead and live cell groups
f = fields(sim_st);
for n=1:numel(f)
    live_st.(f{n}) = sim_st.(f{n})(:,sim_st.isLive);
    dead_st.(f{n}) = sim_st.(f{n})(:,~sim_st.isLive);
end

% keyboard

save('random_exp_sim_st.mat');
load('random_exp_sim_st.mat');

% Calculate ttest
shuffle_inds = randi(tot_trials, [n_resample numel(exp_st.img_name)]);
warning('off','adtest:OutOfRangePHigh');
warning('off','stats:adtest:OutOfRangePHigh');
warning('off','stats:adtest:OutOfRangePLow');
for n = 1:n_resample
    dead_inds = sub2ind([tot_trials n_dead], ...
        shuffle_inds(n,1:n_dead),1:n_dead);
    live_inds = sub2ind([tot_trials n_live], ...
        shuffle_inds(n,n_dead+1:end),1:n_live);
    
    live_st.coloc_score_isnorm(n) = adtest(live_st.coloc_score_rand(live_inds));   
    live_st.coloc_score_ttest_h(n) = ttest(live_st.coloc_score_rand(live_inds),1);
    
    dead_st.coloc_score_isnorm(n) = adtest(dead_st.coloc_score_rand(dead_inds));
    dead_st.coloc_score_ttest_h(n) = ttest(dead_st.coloc_score_rand(dead_inds),1);
end
warning('on','stats:adtest:OutOfRangePHigh');
warning('on','stats:adtest:OutOfRangePLow');
warning('on','adtest:OutOfRangePHigh');


% Coloc score % random
fprintf('Live Cell Coloc Score is Normal in %.2f%% of Cases\n',...
    sum(live_st.coloc_score_isnorm)/n_resample*100)
fprintf('Live Cell Coloc Score is Random in %.2f%% of Cases\n',...
    sum(~live_st.coloc_score_ttest_h)/n_resample*100)
fprintf('Dead Cell Coloc Score is Normal in %.2f%% of Cases\n',...
    sum(dead_st.coloc_score_isnorm)/n_resample*100)
fprintf('Dead Cell Coloc Score is Random in %.2f%% of Cases\n',...
    sum(~dead_st.coloc_score_ttest_h)/n_resample*100)


% Analyze distribution of random variables
live_rv = mean(live_st.coloc_score_rand,2);
live_mean = mean(live_rv);
live_std = std(live_rv);

dead_rv = mean(dead_st.coloc_score_rand,2);
dead_mean = mean(dead_rv); 
dead_std = std(dead_rv);

keyboard
end

