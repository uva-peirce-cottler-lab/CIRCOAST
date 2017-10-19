function [tbl] = ...
    MCM_CirclePacking(img_dim, num_trials, cell_diam, output_path)

CONTAINED_PACKING=1; 
% output_path = [getappdata(0,'proj_path') '/temp_data/CellPackingTest'];
if isempty(dir(output_path)); mkdir(output_path); end

% num_trials=100;  
% img_dim=[512 512];
% cell_diam=9;
cell_rad = (cell_diam-1)/2;

bw_border = true(img_dim);
bw_border(1:cell_rad,:)=0; bw_border(:,1:cell_rad)=0;
bw_border(end-cell_rad+1:end,:)=0; bw_border(:,end-cell_rad+1:end)=0;

% Edge effects are observed at image border, add padding and then run
% simulation image to remove them.
pad = 5*cell_diam;
sim_dim = img_dim+2*pad;

cell_kernel = strel('disk',ceil((cell_diam-1)/2),0);
exclude_kernel = strel('disk',2*ceil((cell_diam-1)/2),0);

% tbl=table();
imwrite(imresize(cell_kernel.Neighborhood,[128 128],'nearest'),...
    [output_path '/' sprintf('kernel_%i_pix.tif',cell_diam)]);
% return

bw_exclude = exclude_kernel.Neighborhood;
% Kernel change in row and columnf rom center
% delta kernal [row,column]
dkr = [round((size(bw_exclude,1)-1)/2) round((size(bw_exclude,1)-1)/2)];
dkc = [round((size(bw_exclude,2)-1)/2) round((size(bw_exclude,2)-1)/2)];
% Kernel center coordinate
kr = dkr(1)+1;
kc = dkc(1)+1;

pack_rat=zeros(num_trials,1);
tot_cells=zeros(num_trials,1);

% keyboard

parfor_progress(50);
parfor k=1:num_trials
    % Reset values
    bw_cell_ctr = false(sim_dim);
    bw_spawn_locs = true(sim_dim);
    img_ind_all = 1:prod(sim_dim);
  
    % Make [row column] matrix coordinates for changing from sub to linear
    % index, done inside parfor because it would be broadcast variable
    r_ind=bsxfun(@times, (1:sim_dim(1))',ones(sim_dim));
    c_ind=r_ind';


    for n=1:numel(bw_spawn_locs)
        % Identify spawn locations in image
        spawn_inds = img_ind_all(bw_spawn_locs);
        if isempty(spawn_inds); break; end
        
        spawn_ind = spawn_inds(randi([1 numel(spawn_inds)],[1 1]));
        % Add all cells to center image
        bw_cell_ctr(spawn_ind)=1;
        
        % Add exclude zone around cell by dilating subimage
%         [r, c] = ind2sub(size(bw_cell_ctr),spawn_ind);
        r = r_ind(spawn_ind);
        c = c_ind(spawn_ind);
        

        % Change in [row,col] coordinates that are allowed in target image
        delt_r = [r-max([1 r-dkr(1)])...
            min([size(bw_spawn_locs,1) r+dkr(2)])-r];
        delt_c = [c-max([1 c-dkc(1)])...
            min([size(bw_spawn_locs,2) c+dkc(2)])-c];
        
        % Fast method of restricting spawn points
        % Since only one point dilated per round, do it manually to save
        % computation time
        bw_spawn_locs([r-delt_r(1):r+delt_r(2)],[c-delt_c(1):c+delt_c(2)]) = ...
            ~(~bw_spawn_locs(r-delt_r(1):r+delt_r(2),c-delt_c(1):c+delt_c(2)) | ...
            bw_exclude(kr-delt_r(1):kr+delt_r(2),kc-delt_c(1):kc+delt_c(2)));
        
%         if mod(n,250)==0; imshow(bw_cell_ctr); end
    end
%     toc
    
    % Excise full image (excluding padded borders with edge effects
    % Excise center image after or before dilation?
    if CONTAINED_PACKING
        % Packing ratio of *fully* contained cells within image borders
        bw_cell_ctr_cont = bw_cell_ctr(pad+1:end-pad,pad+1:end-pad) .* bw_border;
        bw_cells_cont = imdilate(bw_cell_ctr_cont,cell_kernel);
        
    else
        bw_cell_ctr_cont = bw_cell_ctr(pad+1:end-pad,pad+1:end-pad);
        bw_cells = imdilate(bw_cell_ctr_cont,cell_kernel);
        bw_cells_cont= bw_cells(pad+1:end-pad,pad+1:end-pad) 
    end
    
    
    
    % Output variables
    tot_cells(k)=sum(bw_cell_ctr_cont(:));
    pack_rat(k) = sum(bw_cells_cont(:))./numel(bw_cells_cont);
    
    % Write image to file
    pid = feature('getpid');
    imwrite(bw_cell_ctr_cont, [output_path '/bw_ctr_CD' num2str(cell_diam) ...
        '_' num2str(pid) '.tif'],'WriteMode','append');
    
    if mod(k,num_trials/50)==0; parfor_progress; end
    
end
tbl = table(tot_cells, pack_rat,'VariableNames',{'tot_cells','pack_rat'});

% toc
% pack_rat = tot_cells*(pi*cell_diam/2)./prod(img_dim);
% imshow(imdilate(bw_cell_ctr,strel('disk',ceil((cell_diam-1)/2),0)));
