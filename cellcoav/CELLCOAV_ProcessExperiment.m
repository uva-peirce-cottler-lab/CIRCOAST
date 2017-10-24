function CELLCOAV_ProcessExperiment()
% Analyzes cellular colocalization experiment with 2 study groups 
close all

% Flags
FilterUselessImages = 0;
JoinImagesFromGroup = 1;
ClearPreviousData=1;

UseMonteCarloModel=0; 



[csv_name, csv_path] = uigetfile('*.csv','Select CSV of Cell Coloc. Results');
if csv_name==0; return; end
csv_path=csv_path(1:end-1);


if isempty(dir([csv_path '/sim_st.mat'])) || ClearPreviousData
    %% Parse CSV file containing experimental data
    csv_txt = fileread([csv_path '/' csv_name]);
    
    % split text into lines
    csv_lines = regexp(csv_txt, '\r\n|\n\rdb|\n','split');
    
    % Load variables in header of csv file, for example:
    %    BackgroundChannelIndex = 1; group_names={'Live Cells','Dead Cells'};
    elem = @(x) x{1};
    eval(elem(regexp(csv_lines{1},'"(.*)"','tokens','once')));
    
    % Load each column into struct
    csv_fields = cellfun(@(x) regexp(x,',','split'),csv_lines(2:end),'UniformOutput',0);
    nfieldnames = 6;
    csv_fields(cellfun(@(x) numel(x) ~= nfieldnames,csv_fields))=[];
    csv_data = vertcat(csv_fields{:});
    csv_fieldnames = csv_data(1,:);
    for c = 1:size(csv_data,2)
        if isnan(str2double(csv_data{2,c}))
            csv_st.(csv_data{1,c}) = csv_data(2:end,c);
        else csv_st.(csv_data{1,c}) = cellfun(@(x) str2double(x),csv_data(2:end,c));
        end
    end
    for n = 1:numel(csv_st.img_path)
        [~,temp_name] = fileparts(csv_st.img_path{n});
        csv_st.img_name{n} =  [temp_name '.tif'];
    end
    
    % keyboard
    
    % Sort structure by group
    [~,ix] = sort(csv_st.group_id);
    f = fields(csv_st);
    for n=1:numel(f)
        csv_st.(f{n})=csv_st.(f{n})(ix);
    end
    
    % Exclude entries with no cells/events
    ix=csv_st.total==0;
    fprintf('Exluding %i images with no events.\n',sum(ix));
    for n=1:numel(f)
        eval(['csv_st.' f{n} '(ix)=[];'])
    end
    % keyboard
    
    % Identify sample groups so images from same eye are joined
    unq_sample_groups = unique([csv_st.sample_id csv_st.group_id],'rows');
    
    % For each unique sample group, find indiciess for each
    pre_ind = 1:numel(csv_st.sample_id);
    for n = 1:size(unq_sample_groups)
        imgs_ind{n} = pre_ind(unq_sample_groups(n,1) == csv_st.sample_id & ...
            unq_sample_groups(n,2) == csv_st.group_id);
    end
    
    
    if ~JoinImagesFromGroup
        for n = 1:size(unq_sample_groups)
            sample_groups_cell{n}=repmat(unq_sample_groups(n,:),[numel(imgs_ind{n}) 1]);
        end
        unq_sample_groups = vertcat( sample_groups_cell{:});
        imgs_ind=num2cell(horzcat(imgs_ind{:})');
    end
    
    
    % keyboard
    for n = 1:size(unq_sample_groups)
        sim_st.sample_id(n) = unq_sample_groups(n,1);
        sim_st.group_id(n) = unq_sample_groups(n,2);
        
        fprintf('Sample: %.f/%.f, Group %.f/%.f: \n',sim_st.sample_id(n),size(unq_sample_groups,1),...
            sim_st.group_id(n),size(unq_sample_groups,1))

        
        % Load images into cell
        imgs_names='';
        for k = 1:numel(imgs_ind{n})
            img_full_path = [csv_path '\thresh_images\' csv_st.img_name{imgs_ind{n}(k)}];
            img =  imread(img_full_path);% '.tif']);
            if size(img,3)==1; BackgroundChannelIndex=1; end
            imgs_cell{k} = img(:,:,BackgroundChannelIndex)>0;
            
            imgs_names = strcat(imgs_names,csv_st.img_name{imgs_ind{n}(k)},',');
        end
        sim_st.img_vessel_frac(n) = sum(cellfun(@(x) sum(x(:)),imgs_cell))./...
            sum(cellfun(@(x) numel(x(:)),imgs_cell));
        sim_st.imgs_names{n} = imgs_names;
        

        % Calculate stats for hybrid image        
        sim_st.img_names{n} = csv_st.img_name{imgs_ind{n}};
        sim_st.num_imgs(n) = numel(imgs_ind{n});
        sim_st.exp_hits(n) = sum(csv_st.hits(imgs_ind{n}));
        sim_st.exp_misses(n) = sum(csv_st.misses(imgs_ind{n}));
        sim_st.exp_total(n) = sum(csv_st.total(imgs_ind{n}));
        sim_st.exp_frac(n) = sum(csv_st.hits(imgs_ind{n}))/sum(csv_st.total(imgs_ind{n}));
        
        % Get simulated results from Monte Carlo
        if UseMonteCarloModel
            [sim_st.sim_frac_means(n), sim_st.sim_fraction_stds(n), ...
                ~, ~, trials_hit_rate{n}] = ...
                ArcasGui_monteCarloSim_Driver(imgs_cell, cell_diam_um, umppix, ...
                tot_trials, sim_st.exp_total(n), max_trials_per_cycle, handles);
        else
            trials_hit_rate{n}=[];
        end
        
        fract_area = mean(cellfun(@(x) sum(x(:))/numel(x), imgs_cell,'UniformOutput',1));
        fprintf('Observed ICF: %.2f\n MCMRP uICF: %.2f\t%s\n',sim_st.exp_hits(n)/sim_st.exp_total(n),fract_area,'');%imgs_string);
        
        % Calculate mean and std with binomial distribution, and hyp test
        out_st = CELLCOAV_BMRP(imgs_cell, cell_diam_um, umppix,...
            sim_st.exp_total(n),'n_obs_success',sim_st.exp_hits(n));
        f = fields(out_st);
        for k =1:numel(f); eval(['sim_st.' f{k} '(n)=out_st.' f{k} ';']); end
        fprintf('BMRP uICF: %.3f\n', out_st.bmrp_icf_mean);
        
        % Calculate mean and std with binomial distribution, and hyp test
        out_st = CELLCOAV_HMRP(imgs_cell, cell_diam_um, umppix,...
            sim_st.exp_total(n),'n_obs_success', sim_st.exp_hits(n));
        f = fields(out_st);
        for k =1:numel(f); eval(['sim_st.' f{k} '(n)=out_st.' f{k} ';']); end
        fprintf('HMRP uICF: %.3f\n', out_st.hmrp_icf_mean);
        
    end
    

    if FilterUselessImages
        ind = sim_st.exp_total<sim_st.binom_min_ntrials;
        %     ind = sim_st.exp_total >150;
        f = fields(sim_st);
        for n=1:numel(f); sim_st.(f{n})(ind) = []; end
    end
    
    group_ids = unique(sim_st.group_id);
    
    % One sample cellcoav test: does group have enriched coloc.?
    cellcoav_1s_pval = CELLCOAV_1S(sim_st.bmrp_cellcoav_p(sim_st.group_id==group_ids(1)),1);
    legend({'Random',group_names{group_ids(1)}})
    fprintf('[%s] CELLCOAV_1S p:%.4e\n',group_names{group_ids(1)},cellcoav_1s_pval);
    
    cellcoav_1s_pval = CELLCOAV_1S(sim_st.bmrp_cellcoav_p(sim_st.group_id==group_ids(2)),1);
    legend({'Random',group_names{group_ids(2)}})
    fprintf('[%s] CELLCOAV_1S p:%.4e\n',group_names{group_ids(2)},cellcoav_1s_pval);

    
    save([csv_path '/sim_st.mat']);
else
    load([csv_path '/sim_st.mat']);
end


% Plot cell density between groups
inj_cells_p_img = sim_st.exp_total./sim_st.num_imgs;
boxplot(inj_cells_p_img,sim_st.group_id-1)
set(gca,'fontsize',7,'fontname','helvetica')
[h, p]=ttest2(inj_cells_p_img(sim_st.group_id==1),inj_cells_p_img(sim_st.group_id==2));
% mean(inj_cells_p_img(sim_st.group_id==2))./mean(inj_cells_p_img(sim_st.group_id==1))
set(gca,'XTickLabel',group_names,'fontsize',8);
ylabel('Injected Cells/FOV','fontsize',8)
set(gcf,'position', [100 100 160 140])


% EC Cell Density between groups
boxplot(sim_st.img_vessel_frac,sim_st.group_id-1)
set(gca,'fontsize',7,'fontname','helvetica')
[h, p]=ttest2(sim_st.img_vessel_frac(sim_st.group_id==1),sim_st.img_vessel_frac(sim_st.group_id==2));
% mean(sim_st.img_vessel_frac(sim_st.group_id==2))./mean(sim_st.img_vessel_frac(sim_st.group_id==1))
set(gca,'XTickLabel',group_names,'fontsize',7);
ylabel('ECs/FOV','fontsize',8)
set(gcf,'position', [100 100 160 140])


% Two sample cellcoav test, do groups have unique vessel colocalization?
cellcoav_2s_p = CELLCOAV_2S(sim_st.bmrp_cellcoav_p, sim_st.group_id,csv_path,ClearPreviousData);
fprintf('2S CELLCOAV [%s,%s,] p: %.3e\n',...
    group_names{group_ids(1)},group_names{group_ids(2)},cellcoav_2s_p);

 
end

