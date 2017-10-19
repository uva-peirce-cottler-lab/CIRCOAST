function pval = CELLCOAV_2S(ind_var, grp_id, csv_path,ClearPreviousData)
% Generate permuted sample IDs and then save to disk
n_rounds = 1e7;

ugrp_id = unique(grp_id); 



% Compute test statistic on groups
[p_val_actual] = ranksum(ind_var(grp_id==ugrp_id(1)),ind_var(grp_id==ugrp_id(2)));
n_samples = numel(ind_var);

% A large amount of permutations (each a permuted index of N samples) is
% needed for the analysis, If there is a permutation already cached on disk
% with the same number of samples, load instead of recalculate.
permutation_ind_path = [csv_path '/permutation_temp_data.mat'];

% Load previous data for permutations if exist
if ~isempty(dir(permutation_ind_path))
    st = load(permutation_ind_path);
else st=[];
end

% Make sure that n_samples and n_rounds is right.
if isempty(st) || ClearPreviousData ||...
        (size(st.sim_sample_ids,1)~=n_rounds && size(st.sim_sample_ids,2)~=n_samples)
    sim_sample_ids = zeros([n_rounds,n_samples],'uint8');
    
    hw = waitbar(0,'Permuting IDs...');
    for n = 1:n_rounds
        sim_sample_ids(n,:)=grp_id(randperm(n_samples));
        if mod(n,n_rounds/100)==0; waitbar(n/n_rounds,hw); end
    end; close(hw)
    save(permutation_ind_path,'sim_sample_ids','-v7.3');
    
end
if ~exist('sim_sample_ids','var')
    load(permutation_ind_path);
end

%% Pre-sort data into two groups
grp_var1 = zeros(n_rounds,sum(grp_id==ugrp_id(1)));
grp_var2 = zeros(n_rounds,sum(grp_id==ugrp_id(2)));
hw = waitbar(0,'Splitting group data...');
for n=1:n_rounds
   grp_var1(n,:)= ind_var(sim_sample_ids(n,:)==ugrp_id(1));
   grp_var2(n,:)= ind_var(sim_sample_ids(n,:)==ugrp_id(2));
    if mod(n,n_rounds/100)==0; waitbar(n/n_rounds,hw); end
end; close(hw);

% Randomly permute study group number assignment, create distribution of test statistic
%Cache disk with premutation IDs
permutation_pval_path = [csv_path '/temp_data/permutation_pvalues.mat'];
mkdir([csv_path '/temp_data/']);
if isempty(dir((permutation_pval_path))) || ClearPreviousData
    
    % Stare pvalues in var
    p_val_sim = NaN(1,n_rounds);
    
    % Get number of logical cores
    num_log_cores = str2double(regexp(evalc('feature(''numcores'')'),...
        '(\d*).logical','tokens','once'));
    
    %Start parallel processing
    delete(gcp)
    c = parcluster;
    c.NumWorkers = max([num_log_cores-1 3]);
    saveProfile(c);
    parpool(7);
    
    fprintf('%s\n',datetime('now'));
    parfor_progress(100);
    parfor n=1:n_rounds
         [p_val_sim(n)] = ranksum(grp_var1(n,:), grp_var2(n,:));
        if mod(n,n_rounds/100)==0; parfor_progress; end
    end
    parfor_progress(0);
    fprintf('%s\n',datetime('now'));
    

    save(permutation_pval_path,'p_val_sim','-v7.3');
else
    %     keyboard
    load(permutation_pval_path);
    
end
% toc

% Grab mean and std of distribution of p values
stats.mean_grp_pval = mean(p_val_sim);
stats.std_grp_pval= std(p_val_sim);

% Sort distribution and find percentile where actual p value falls
sorted_p_val_sim = sort(p_val_sim);

elem = @(x,k) x(k);
% See where test statistic with real group #
bv = sorted_p_val_sim>(p_val_actual);%-0.5);
ind = elem(elem(1:n_rounds,bv),1);
prct_permute_1t_p_mean = ind/n_rounds;
pval=prct_permute_1t_p_mean;


%Plot Results
figure;
hist_data = histogram(sorted_p_val_sim,25,'Normalization','probability');
% beautifyAxis(gca);
hold on
y=ylim; x=[0 1];
plot([prct_permute_1t_p_mean prct_permute_1t_p_mean],...
    [0 y(2)],'LineWidth',2.5,'Color','R')
% beautifyAxis(gca);
legend({'Random','Observed'})
xlabel('Permuted CELLCOAV p Values')
ylabel('Probability')
% title(sprintf(['Permutation Test of Study Groups 1 and 2 \n with Wilsox Sum Rank Test of Binomial p Values']))
% legend({'Permuted P Value','Actual P Value'})
set(findall(gcf,'-property','FontSize'),'FontSize',7);
set(findall(gcf,'-property','FontName'),'FontName','Helvetica');
set(gcf,'position', [100 100 200 150])
axis([x y]);
hold off

fprintf('GroupPerm P value: %f\n', pval);

keyboard
