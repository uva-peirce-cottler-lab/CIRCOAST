function [] = SimExp_ImageJoining()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
EXPORT_IMAGES=0;
ClearPreviousData=1;
%Output Directory
output_path = [getappdata(0,'proj_path') '/temp_data/ImageJoinTest'];
mkdir(output_path);
% delete([output_path '/*.*']);

% Simulation Parameters
N_tot=500;
img_dim = [512 512];
umppix = 424.5/512;
cell_diam_um = 8;
cell_diam_pix = cell_diam_um/umppix;


% Data table
% image_Name, group_id, sample_id, rep_id, vld, tot_cells, obs_col_cells, obs_icf, bmrp_icf
tbl=table;
tbl.img_name = cell(N_tot,1);
tbl.vld = zeros(N_tot,1);
tbl.tot_cells = round(normrnd(ones(N_tot,1).*80,5));
tbl.whole_img_obs_icf = zeros(N_tot,1);
tbl.whole_img_mean_icf = zeros(N_tot,1);
% tbl.whole_img_cellcoav_p = zeros(N_tot,1);
tbl.joined_img_mean_icf = zeros(N_tot,1);
% tbl.joined_img_cellcoav_p = zeros(N_tot,1);


% Requested vld is not exactly the same as actual
% Save time by regressing 4 additional times after each vessel network
% instead of producing from scratch
sub_vect = bsxfun(@times, rand(5,N_tot/5), [0; 5; 10; 15; 20]);
req_vld = normrnd(vertcat(ones(N_tot,1).*35),2)-sub_vect(:);


% Create image dataset
hw = waitbar(0,'Processing Images');
bw_skel = false(img_dim);
k=1;
for N=1:N_tot
    tbl.img_name{k} = sprintf('I%i.tif',N);
    
    if isempty(dir([output_path '/' tbl.img_name{k}]))
        % If bw_skel has less VLD than target, make new vessel network
        if sum(bw_skel(:))*(umppix/1000)/(img_dim(1)*umppix/1000).^2 < req_vld(N)
            
            % Create Initial over dense vessel network
            [~, init_bw_skel, ~] = VesselGen_GenerateHoneycombNetwork(img_dim, 3,umppix,...
                3*6, 'VesselLengthDensityLimits', [50 0],'FigureVisible','on');
            % Identify central vessel to protect
            bw_protect = VesselGen_SelectProtectedSegment(init_bw_skel,'FigureVisible','on');
        end
        
        % Iteratively remove vessel segments to a range of vessel densities
        [bw_vessel, bw_skel, stat_st] =...
            VesselGen_RegressNetwork(init_bw_skel, 3,umppix, ...
            'VesselLengthDensityTarget',req_vld(N),...
            'BW_Protected',bw_protect,'FigureVisible','on');
        tbl.vld(N) = stat_st.vld_mmpmm2;
        
        
        % Distribute cells in image
        [~, ~, ~,comp_img,tbl.whole_img_obs_icf(N)] = ...
            CELLCOAV_MCMRP(bw_vessel, 1, ...
            tbl.tot_cells(N), cell_diam_pix);
        
        % Split image into 4 parts
        sub_img = quarter_split_img(comp_img);
        sub_vasc_img = cellfun(@(x) x(:,:,2),sub_img,'UniformOutput',0);
        
        
        % Export Images
        tbl.img_name{N} =  sprintf('I%i.tif',N);
        imwrite(comp_img,[output_path '/' tbl.img_name{N}]);
        for n=1:4
            imwrite(sub_img{n},[output_path '/' sprintf('I%i_%i.tif',N,n)]);
        end
    else
        % Load Images
        comp_img = imread([output_path '/' tbl.img_name{k}]);
        for n=1:4
            sub_img{n} = imread([output_path '/' sprintf('I%i_%i.tif',N,n)]);
        end
        sub_vasc_img = cellfun(@(x) x(:,:,2),sub_img,'UniformOutput',0);
        
        tbl.vld(N) = skel_2_vld_mmpmm2(bwmorph(comp_img(:,:,2),'Thin',Inf), umppix);
        tbl.tot_cells(N) = sum(sum(comp_img(:,:,3)>0));
        bw_vessel = comp_img(:,:,2);
        
        %         keyboard
    end
    
    tbl.coloc_cells(N) = sum(sum(comp_img(:,:,3)==2));
    
    % Calculate ICF for entire image
    binom_st = CELLCOAV_BMRP(bw_vessel, ...
        cell_diam_um, umppix,tbl.tot_cells(N),'n_obs_success',...
        tbl.coloc_cells(N));
    % Get p value form whole image
    tbl.whole_img_mean_icf(N) = binom_st.bmrp_icf_mean;
    tbl.whole_img_cellcoav_p(N) = binom_st.bmrp_cellcoav_p;
    % Calculate ICF for joined images
    binom_st = CELLCOAV_BMRP(sub_vasc_img, ...
        cell_diam_um, umppix,tbl.tot_cells(N),'n_obs_success',...
        tbl.coloc_cells(N));
    % Get p value form joined images
    tbl.joined_img_mean_icf(N) = binom_st.bmrp_icf_mean;
    tbl.joined_img_cellcoav_p(N) = binom_st.bmrp_cellcoav_p;
    
    
    % Clear border of image and recalculate
    bord_comp_img = comp_img;
    bord_comp_img(round(img_dim(1)/2)-ceil(cell_diam_pix/2): ...
        round(img_dim(1)/2)+ceil(cell_diam_pix/2),:,2)=0;
    bord_comp_img(:,round(img_dim(2)/2)-ceil(cell_diam_pix/2): ...
        round(img_dim(2)/2)+ceil(cell_diam_pix/2),2)=0;
    % Split image into 4 parts
    bord_sub_img = quarter_split_img(bord_comp_img);
    bord_sub_vasc_img = cellfun(@(x) x(:,:,2),bord_sub_img,'UniformOutput',0);
    
    if EXPORT_IMAGES
        % export Border Images
        imwrite(bord_comp_img,[output_path '/bord_' tbl.img_name{N}]);
        for n=1:4
            imwrite(bord_sub_img{n},[output_path '/bord_' sprintf('I%i_%i.tif',N,n)]);
        end
    end
    
    % Calculate ICF for entire image
    binom_st = CELLCOAV_BMRP(bord_comp_img(:,:,2), ...
        cell_diam_um, umppix,tbl.tot_cells(N),'n_obs_success',...
        tbl.coloc_cells(N));
    % Get p value form whole image
    tbl.bord_whole_img_mean_icf(N) = binom_st.bmrp_icf_mean;
    tbl.bord_whole_img_cellcoav_p(N) = binom_st.bmrp_cellcoav_p;
    % Calculate ICF for joined images
    binom_st = CELLCOAV_BMRP(bord_sub_vasc_img, ...
        cell_diam_um, umppix,tbl.tot_cells(N),'n_obs_success',...
        tbl.coloc_cells(N));
    % Get p value form joined images
    tbl.bord_joined_img_mean_icf(N) = binom_st.bmrp_icf_mean;
    tbl.bord_joined_img_cellcoav_p(N) = binom_st.bmrp_cellcoav_p;
    
    
    %             keyboard
    waitbar(k/(N_tot),hw);
    k=k+1;
    
end

% tbl.coloc_cells = tbl.whole_img_obs_icf .* tbl.tot_cells;

close(hw);

% Export table
writetable(tbl,[output_path '/arcas_results.csv']);
save([output_path '/arcas_data.m']);

keyboard
% clear all
output_path = [getappdata(0,'proj_path') '/temp_data/ImageJoinTest'];
load([output_path '/arcas_data.m']);
tbl=readtable([output_path '/arcas_results.csv']);
%% Pval comparison
figure;
[h,p] = ttest(tbl.whole_img_cellcoav_p, tbl.joined_img_cellcoav_p);
x=(tbl.whole_img_cellcoav_p + tbl.joined_img_cellcoav_p)./2;
y=tbl.whole_img_cellcoav_p - tbl.joined_img_cellcoav_p;
fprintf('TTest for Whole vs. Joined CELLCOAV p: %.4e\n',p);
plot(x,y,'ro', 'MarkerSize', 3);
[muhat,sigmahat,muci,sigmaci] = normfit(tbl.whole_img_mean_icf - tbl.joined_img_mean_icf);
hold on;
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y) mean(y)+1.96*std(y)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b');
hold off
xlabel('\mu CELLCOAV p of Whole and Split Image');
ylabel('Whole \mu CELLCOAV p - Split \mu CELLCOAV p')
beautifyAxis(gcf);
ya_pval=ylim;
set(gcf,'position',[300 300 300 250])

figure;
[h,p] = ttest(tbl.bord_whole_img_cellcoav_p, tbl.bord_joined_img_cellcoav_p);
x=(tbl.bord_whole_img_cellcoav_p + tbl.bord_joined_img_cellcoav_p)./2;
y=tbl.bord_whole_img_cellcoav_p - tbl.bord_joined_img_cellcoav_p;
fprintf('TTest for Border Whole vs. Joined CELLCOAV p: %.4e\n',p);
plot(x,y,'ro', 'MarkerSize', 3);
[muhat,sigmahat,muci,sigmaci] = normfit(tbl.whole_img_mean_icf - tbl.joined_img_mean_icf);
hold on;
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y) mean(y)+1.96*std(y)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b');
hold off
xlabel('\mu CELLCOAV p of Whole and Split Image');
ylabel('Whole \mu CELLCOAV p - Split \mu CELLCOAV p')
beautifyAxis(gcf);
axis([xlim, ya_pval-mean(ya_pval)]);
set(gcf,'position',[300 300 300 250])

%% CDVF comparison
figure;
[h,p] = ttest(tbl.whole_img_mean_icf, tbl.joined_img_mean_icf);
x=(tbl.whole_img_mean_icf + tbl.joined_img_mean_icf)./2;
y=tbl.whole_img_mean_icf - tbl.joined_img_mean_icf;
fprintf('TTest for Whole vs. Joined ICF: %.4e\n',p);
plot(x,y,'ro', 'MarkerSize', 3);
[muhat,sigmahat,muci,sigmaci] = normfit(tbl.whole_img_mean_icf - tbl.joined_img_mean_icf);
hold on;
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y) mean(y)+1.96*std(y)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b');
hold off
xlabel('\mu CDVG of Whole and Split Image');
ylabel('Whole \mu CDVF - Split \mu CDVF');
beautifyAxis(gcf);
ya_CDVF=ylim;
set(gcf,'position',[300 300 300 250])

figure;
[h,p] = ttest(tbl.bord_whole_img_mean_icf, tbl.bord_joined_img_mean_icf);
x=(tbl.bord_whole_img_mean_icf + tbl.bord_joined_img_mean_icf)./2;
y=tbl.bord_whole_img_mean_icf - tbl.bord_joined_img_mean_icf;
fprintf('TTest for Border Whole vs. Joined ICF: %.4e\n',p);
plot(x,y ,'ro', 'MarkerSize', 3);
hold on;
plot([xlim' xlim'], repmat([mean(y)-1.96*std(y) mean(y)+1.96*std(y)],[2 1]),'g');
plot(xlim, [mean(y) mean(y)],'b');
hold off
xlabel('\mu CDVF of Whole and Split Image');
ylabel('Whole \mu CDVF - Split \mu CDVF');
beautifyAxis(gcf);
axis([xlim, ya_CDVF-mean(ya_CDVF)]);
set(gcf,'position',[300 300 300 250])

keyboard
end


function sub_img = quarter_split_img(comp_img)
img_dim = size(comp_img);
% Get p value from images split 4 times
sub_img{1}= comp_img(1:floor(img_dim(1)/2),1:floor(img_dim(2)/2),:);
sub_img{2}= comp_img(1:floor(img_dim(1)/2),floor(img_dim(2)/2)+1:img_dim(2),:);
sub_img{3}= comp_img(floor(img_dim(1)/2)+1:img_dim(1),1:floor(img_dim(2)/2),:);
sub_img{4}= comp_img(floor(img_dim(1)/2)+1:img_dim(1),...
    floor(img_dim(2)/2)+1:img_dim(2),:);
end
