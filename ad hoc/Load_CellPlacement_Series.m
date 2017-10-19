function zbw_ctrs  = ...
    Load_CellPlacement_Series(rexp, img_path,SPLIT)

% keyboard
item_list = dir([img_path '/*']);
item_cell = {item_list(:).name}';

ind =1:numel(item_cell);
bv_matches = cellfun(@(x) ~isempty(regexp(x, rexp,'once')),item_cell);
bv_ind = ind(bv_matches);

img_cell = {numel(bv_ind),1};

for n=1:numel(bv_ind)
   
    img_info=imfinfo([img_path '/' item_cell{bv_ind(n)}]);
    img_cell{n}=false(img_info(1).Height,img_info(1).Width,length(img_info));
    for z=2:2:numel(img_info)
        img_cell{n}(:,:,z)=imread([img_path '/'  item_cell{bv_ind(n)}],'Index',z,...
            'Info',img_info);
    end
    
end

zbw_ctrs = cat(3,img_cell{:});


% zbw_cells=zimg(:,:,1:2:end);
% zbw_ctrs = zimg(:,:,2:2:end);












