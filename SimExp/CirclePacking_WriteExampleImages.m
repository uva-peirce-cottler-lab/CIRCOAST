function CirclePacking_WriteExampleImages(img_path)


bw_cell_starr = dir([img_path '/bw_cell_*.tif']);

bw_cell_names = {bw_cell_starr(:).name};

cell_diams = cellfun(@(x) str2double(x), ...
    regexp(bw_cell_names, 'bw_cell_CD(\d*)','once','tokens'));

unq_cell_diams = unique(cell_diams);

for n=1:numel(unq_cell_diams)
   img_name = bw_cell_names{find(unq_cell_diams(n)==cell_diams,1)};
   bw = imread([img_path '/' img_name]);
    
   bw_cell = imdilate(bw, strel('disk',(unq_cell_diams(n)-1)/2,0));
   
   imwrite(bw_cell, [img_path '/' sprintf('bw_cells_CD%i.tif',unq_cell_diams(n))]);
end