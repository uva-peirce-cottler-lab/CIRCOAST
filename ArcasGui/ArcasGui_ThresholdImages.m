
function ArcasGui_ThresholdImages(imgs_path,thesh_chan_index)
% imgs_path
base_path = fileparts(imgs_path);


if isempty(dir([base_path '/thresh_images/'])); mkdir([base_path '/thresh_images/']); end
items = dir([imgs_path '\*.tif']);


% Eliminate imgs with _bv prefix
imgs_name={items(:).name};
% bw_img_names = imgs_name(cellfun(@(x) ~isempty(regexp(x,'.tif_bw.tif', 'once')),imgs_name));
% 
% for n = 1: numel(bw_img_names)
%     delete([curr_path '/' bw_img_names{n}]);
% end

imgs_name(cellfun(@(x) ~isempty(regexp(x,'.tif_bw.tif', 'once')),imgs_name))=[];

fprintf('Processing %.f: %s\n', numel(imgs_name), imgs_path);


for n = 1:numel(imgs_name)
    fprintf('%s\n',imgs_name{n});
   rgb_img=imread([imgs_path '/' imgs_name{n}]);
   imshow(imresize(rgb_img,[512 512])); pause(.1);
   if size(rgb_img,3)==1; 
       img = rgb_img(:,:,1); 
   else
       img = rgb_img(:,:,thesh_chan_index);  
   end 
   
   bimg = medfilt2(img,[200 200],'symmetric');
   
   fimg = img-bimg> 40;
   
%    bfimg = bwareaopen(fimg, 50);
   bfimg= img>0;
   imshow(imresize(bfimg,[512 512]));
   
  ind=1:3; ind(thesh_chan_index)=[];
  rgb_img(:,:,ind)=0;
%   keyboard
    imwrite(rgb_img,[base_path '/thresh_images/' regexprep(imgs_path(numel(base_path)+2:end),'\\|/','_') ...
       '_' items(n).name]);
   imwrite(bfimg,[base_path '/thresh_images/' regexprep(imgs_path(numel(base_path)+2:end),'\\|/','_') ...
       '_' items(n).name '_bw.tif']);
end


% keyboard
items = dir([imgs_path '\*']);
items(~cell2mat({items(:).isdir}))=[];
items(1:2)=[];
items(arrayfun(@(x) ~isempty(regexp(x.name,'thresh_images', 'once')),items))=[];

for n = 1:numel(items)
    ArcasGui_ThresholdImages([imgs_path '/' items(n).name],base_path,thesh_chan_index);
end


end