clear all
close all

CELL_AREA_FRACT=1;
% dilCellFcn= @(x) ArcasGui_PerccentCellOverlap(imgs_cell,cell_diam_um,umppix)

% Calculate simulated mean colofcalization fraction with a monte carlo and
% binomial model of random placement and plot results
% Parameters
cell_diam_um = 10;
umppix = 424.5/512;


img_dim = [512 512];
tot_cells = 20;
tot_trials=1e6;


fcn_srt_rand_vessel_density = @(x) sort(randi(23,[1 x])+7,'descend');
fcn_rand_vessel_rad_pix = @(x) randi(20,[1 x])+1;
fcn_rand_cell_diam_um = @(x)  randi(33,[1 x])+2;
fcn_rand_cell_num = @(x)  randi(39,[1 x])+1;

proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/VF_vs_VDF'];
if isempty(dir(out_path)); mkdir(out_path); end


delete([out_path '/*.tif']);

% 10 image series, 20 cells
% Each with
%    5 different vessel desnities
%    5 different vessel_rads
%    5 Different cell sizes

% Qantify mean coloc frac
% VLD ,VF, Dil VF, Cell Diam, Vesel_Diam


% permSampleOrdered = @(x,n) x(sort(randperm(numel(x),n),'ascend'));
% permSample = @(x,n) x(randperm(numel(x),n));

pred_st=struct();
n_p_img = 50;

multiWaitbar('Images #', 0 );
k=1;
for r = 1:10
    
    vect_vessel_len_density = fcn_srt_rand_vessel_density(5);
    
    %     restrict_rad = range_vessel_rad_pix(2)*6;
    [~, bw_skel, stat_st] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
        3*6, 'VesselLengthDensityLimits', [50 0]);
    
    % Identify central vessel to protect
    bw_protect = VesselGen_SelectProtectedSegment(bw_skel);
    
    % Iteratively remove vessel segments to a range of vessel densities
    [~, skel_cell{1}, stat_st] =...
        VesselGen_RegressNetwork(bw_skel, 3,umppix, ...
        'VesselLengthDensityTarget',vect_vessel_len_density(1),'BW_Protected',bw_protect);
    for c=2:5
        [~, skel_cell{c}, stat_st] =...
            VesselGen_RegressNetwork(skel_cell{c-1}, 3,umppix, ...
            'VesselLengthDensityTarget',vect_vessel_len_density(c),'BW_Protected',bw_protect);
        
    end
    
    
    
    
    for c=1:5
        vect_vessel_rad_pix = fcn_rand_vessel_rad_pix(n_p_img);
        vect_cell_diam_um = fcn_rand_cell_diam_um(n_p_img);
        vect_cell_num = fcn_rand_cell_num(n_p_img);
%         vecto_tot_cells = fcn_rand_cell_num(n_p_img);
        for nn=1:n_p_img
            
            cell_diam_um = vect_cell_diam_um(nn);
            cell_num = vect_cell_num(nn);
            vessel_rad_pix=vect_vessel_rad_pix(nn);
            
            % Grab image, dilate, reskeletonize
            bw_vessel = bwdist(skel_cell{c}) <=vessel_rad_pix;
            bw_skel = bwmorph(bw_vessel,'skel',Inf);
            
            % Get stats
            pred_st.cell_diam_um(k) = cell_diam_um;
            pred_st.cell_num(k) = cell_num;
            pred_st.vessel_rad_um(k) = vessel_rad_pix.*umppix;
            pred_st.vessel_frac(k) = sum(bw_vessel(:))./prod(img_dim);
            pred_st.vld_mmpmm2(k) = ImageQuant_Pix2VesselLengthDensity(bw_skel,umppix);
            pred_st.dil_vessel_frac(k)=CELLCOAV_CellDilateVessFrac(bw_vessel,...
                cell_diam_um,umppix);
            
            
            % Compute Response variable
            [resp_st.mcm_colfrac_means(k), resp_st.mcm_colfract_stds(k), ~, comp_img, ~] = ...
                ArcasGui_monteCarloSim_Driver(bw_vessel, cell_diam_um, umppix, ...
                tot_trials, cell_num, 10000);
            
            
            % Get binomial stats
            binom_st = CELLCOAV_BMRP(bw_vessel, ...
                cell_diam_um, umppix,cell_num);
            resp_st.bmd_colfrac_means(k) = binom_st.binom_frac_mean;
            resp_st.bmd_colfrac_stds(k) = binom_st.binom_frac_std;
            
            
            % Randomly generate n binomial random fractions to match MCM
            brf = binornd(cell_num,binom_st.binom_frac_mean, [1 tot_trials])/tot_cells;
            resp_st.bms_colfrac_means(k) = mean(brf);
            resp_st.bms_colfrac_stds(k) = std(brf);
            
            
    
            
            % Output example image
            cell_ind = randi(prod(img_dim),[1 tot_cells]);
            bw_cell=false(img_dim);
            bw_cell(cell_ind)=1;
            
            bw_example = cat(3,bwdist(bw_cell)<=pred_st.cell_diam_um(k),bw_vessel,false(img_dim));
            imwrite(im2uint8(bw_example), [out_path '/' sprintf('%i',k) '.tif'])
            
            multiWaitbar('Images #', k/(10*5*50) );
            k=k+1;
        end
    end
end

save([out_path '/coloc_fraction_predictors.mat'])


% Run here to load data from last run
proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/ParSweep/VF_vs_VDF'];
load([out_path '/coloc_fraction_predictors.mat'])

keyboard

% Convert vessel pix to um
% pred_st.vessel_rad_um=pred_st.vessel_rad_pix
% pred_st=rmfield(pred_st,'vessel_rad_pix')

% Get predictor values
fe = fields(pred_st);
pred_mat= cell2mat(struct2cell(pred_st))';


str= char(fe')';


%Write to csv file
f=fopen([out_path '/coloc_fraction_predictors.csv'],'w');
fprintf(f,'%s,response,\n', regexprep(str(:)','\s*',','));
for n=1:size(pred_mat,1)
    %     sprrintf(
    fprintf(f,'%s,%.4f,\n', regexprep(num2str(pred_mat(n,:)),'\s*',','),...
        resp_st.mcm_colfrac_means(n));
end
fclose(f);

plot(pred_mat,bsxfun(@times,resp_st.bmd_colfrac_means',ones(size(pred_mat))),'x')

pred_st
xtext={'Cell Diam (um)','Cell Number','Vessel Radius (um)','Vessel Fraction',...
    'Vessel Len. Dens. (mm/mm2)','Cell Dil. Vessel Frac'};

close all
for n=1:numel(fe)
    figure;
    plot(pred_mat(:,n), resp_st.bmd_colfrac_means,'r.')
    xlabel(xtext{n});
    ylabel('Mean Cell Coloc. Fraction')
    legend('MCM-RP')
    [RHO(n),PVAL(n)] = corr(pred_mat(:,n),resp_st.mcm_colfrac_means');
    beautifyAxis(gca)
    %    legend('BMRP')
    %     set(gcf,'Position', [100 100 280 260])
    %    xh = get(gca,'XLabel');
    %    set(xh, 'Units', 'Normalized');
    %    set(xh,'Position',get(xh,'Position') + [0 .02 0]);
    set(gcf,'Position', [100 100 250 280])
    
    saveas(gcf,[out_path '/' regexprep(xtext{n},'\/','_') '.fig'])
    
    mdl = fitlm(pred_mat(:,n), resp_st.bmd_colfrac_means);
    
    mRsq(n) = dl.Rsquared.Ordinary
    
    pause();
    close(gcf)
end

close all
% b = regress(resp_st.bmd_colfrac_means',pred_mat);


%% Lasso Regression
[B,FitInfo] = lasso(pred_mat,resp_st.mcm_colfrac_means');
lassoPlot(B,FitInfo);

%% Step Wise REgression
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(pred_mat,resp_st.mcm_colfrac_means');


%% Linear Model Fit
mdl = fitlm(pred_mat,resp_st.mcm_colfrac_means','Linear','RobustOpts','on',...
    'Intercept',false,'PredictorVars',xtext);
anova(mdl,'summary')
stats = anova(mdl)
[c,~,~,gnames] = multcompare(stats);
coefCI(mdl)
[p,F,d] = coefTest(mdl);
plot(mdl)
beautifyAxis(gca)
set(gcf,'Position', [100 100 320 260])


%% Model Comparison
[h,p] = ttest(resp_st.mcm_colfrac_means,resp_st.bmd_colfrac_means)

boxplot(pred_mat,'orientation','horizontal','labels',xtext)
C = corr(pred_mat,pred_mat);

[wcoeff,score,latent,tsquared,explained] = pca(pred_mat,...
'VariableWeights','variance');


% [coeff,score,latent,tsquared,explained,mu] = pca(pred_mat);

pointsize = 10;
scatter(score(:,1),score(:,2), pointsize, resp_st.bmd_colfrac_means-resp_st.mcm_colfrac_means);
colormap hsv
colorbar

xlabel('Principle Comp. 1')
ylabel('Principle Comp. 2')
beautifyAxis(gca)



