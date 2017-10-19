clear all
close all


proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/CellPlacementModel_Comparison'];
if isempty(dir(out_path)); mkdir(out_path); end


CELL_AREA_FRACT=1;
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

pred_tbl=table();
pred_tbl.img_index = (1:2500)';
pred_tbl.cell_diam_um=zeros(2500,1);
pred_tbl.cell_num = zeros(2500,1);
pred_tbl.vessel_rad_um = zeros(2500,1);
pred_tbl.vessel_frac = zeros(2500,1);
pred_tbl.vld_mmpmm2 = zeros(2500,1);
pred_tbl.dil_vessel_frac = zeros(2500,1);

resp_tbl=table();
resp_tbl.img_index = (1:2500)';
resp_tbl.mcmrp_icf_mean=zeros(2500,1);
resp_tbl.mcmrp_icf_std=zeros(2500,1);
resp_tbl.bmrp_icf_mean=zeros(2500,1);
resp_tbl.bmrp_icf_std=zeros(2500,1);
resp_tbl.hmrp_icf_mean=zeros(2500,1);
resp_tbl.hmrp_icf_std=zeros(2500,1);

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
            pred_tbl.cell_diam_um(k) = cell_diam_um;
            pred_tbl.cell_num(k) = cell_num;
            pred_tbl.vessel_rad_um(k) = vessel_rad_pix.*umppix;
            pred_tbl.vessel_frac(k) = sum(bw_vessel(:))./prod(img_dim);
            pred_tbl.vld_mmpmm2(k) = ImageQuant_Pix2VesselLengthDensity(bw_skel,umppix);
            pred_tbl.dil_vessel_frac(k)=CELLCOAV_CellDilateVessFrac(bw_vessel,...
                cell_diam_um,umppix);
            
            
            % Compute Response variable
            [resp_tbl.mcmrp_icf_mean(k), resp_tbl.mcmrp_icf_std(k), ~, comp_img, ~] = ...
                ArcasGui_monteCarloSim_Driver(bw_vessel, cell_diam_um, umppix, ...
                tot_trials, pred_tbl.cell_num(k), 10000);
            
            
            % Get binomial stats
            binom_st = CELLCOAV_BMRP(bw_vessel, ...
                cell_diam_um, umppix,pred_tbl.cell_num(k));
            resp_tbl.bmrp_icf_mean(k) = binom_st.binom_frac_mean;
            resp_tbl.bmrp_icf_std(k) = binom_st.binom_frac_std;
             
            % Get hypergeometric stats
            hypg_st = CELLCOAV_HMRP(bw_vessel, ...
                cell_diam_um, umppix, pred_tbl.cell_num(k));
            resp_tbl.hmrp_icf_mean(k) =hypg_st.hmrp_icf_mean;
            resp_tbl.hmrp_icf_std(k) = hypg_st.hmrp_icf_std;
            
            
            % Output example image
            cell_ind = randi(prod(img_dim),[1 tot_cells]);
            bw_cell=false(img_dim);
            bw_cell(cell_ind)=1;
            
            bw_example = cat(3,bwdist(bw_cell)<=pred_tbl.cell_diam_um(k),bw_vessel,false(img_dim));
            imwrite(im2uint8(bw_example), [out_path '/' sprintf('%i',k) '.tif'])
            
            multiWaitbar('Images #', k/(10*5*50) );
%             keyboard
            k=k+1;
        end
    end
end

save([out_path '/coloc_fraction_predictors.mat'])


% Run here to load data from last run
proj_path = getappdata(0, 'proj_path');
out_path = [proj_path '/temp_data/CellPlacementModel_Comparison'];
load([out_path '/coloc_fraction_predictors.mat'])


% Recrunch hypergeomtric
% Load each image
% recrunch and renter info in table.



%Write pred and response to csv file
writetable(pred_tbl,[out_path  '/predictor_vars.csv']);
writetable(resp_tbl,[out_path  '/presponse_vars.csv']);


% Labels of each predictor
pred_label={'Cell Diam (um)','Cell Number','Vessel Radius (um)','Vessel Fraction',...
    'Vess. Len. Dens. (mm/mm2)','Cell Dil. Vessel Frac'};

%% MCMRP vs Preditors

close all
fe = pred_tbl.Properties.VariableNames;
fe(1)=[];
% Plot MCMRP mean ICF over predictor variables
fprintf('Pearson Test of MCMRP vs Predictors\n')
for n=1:numel(fe)
    figure;
    plot(pred_tbl.(fe{n}), resp_tbl.mcmrp_icf_mean,'r.','MarkerSize',4)
    if strcmp(fe{n},'dil_vessel_frac'); hold on
       plot([0 1],[0 1],'k');hold off 
    end
    xlabel(pred_label{n});
    ylabel('Mean ICF')

    [R,P,RLO,RUP]= corrcoef(pred_tbl.(fe{n}), resp_tbl.mcmrp_icf_mean, 'alpha', 0.05);
    fprintf('%s (r=%.2f [%.2f %.2f] %s)\n',pred_label{n},...
        R(1,2),RLO(1,2),RUP(1,2),printf_pval(P(1,2)));
    beautifyAxis(gca);   
    set(gcf,'Position', [100 100 155 175]);
    set(findall(gcf,'-property','FontSize'),'FontSize',8);
    set(findall(gcf,'-property','FontName'),'FontName','Helvetica');
    set(gca,'XTickLabel',get(gca,'XTickLabel'),'fontsize',7);
    saveas(gcf,[out_path '/' regexprep(pred_label{n},'\/','_') '.fig']);
    pause();
    close(gcf)
end


%% Import data from spreadsheet instead of data folder
% tbl = readtable([out_path '/ICF_ModelComp.csv']);
% pred_tbl(:,2:7)=tbl(:,1:6);
% resp_tbl.mcmrp_icf_mean=tbl.mean_mcmrp_icf;
% resp_tbl.bmrp_icf_mean=tbl.mean_bmrp_icf;

%% MCMRP to PREDICTORS: Linear Model Fit
% Calculate z scores for predictor to make comparison valid
% Zscore of predictors: (x-mean(x))/std(x)
mdl = fitlm(zscore_col(table2array(pred_tbl(:,2:end))),zscore_col(resp_tbl.mcmrp_icf_mean),'Linear',...
    'Intercept',true,'PredictorVars',pred_label);
mdl.Coefficients.SE
% http://stats.stackexchange.com/questions/79399/calculate-variance-explained-by-each-predictor-in-multiple-regression-using-r
MCMRP_Annova = anova(mdl);




SimExp_CellPlacementModel_Comparison_Plots(pred_tbl,resp_tbl,'MCMRP','mcmrp_icf_mean', ...
    'BMRP','bmrp_icf_mean', pred_label,out_path);


SimExp_CellPlacementModel_Comparison_Plots(pred_tbl,resp_tbl,'BMRP','bmrp_icf_mean', ...
    'HMRP','hmrp_icf_mean', pred_label,out_path);


