function [bw_imgs, bw_rois] = ArcasGui_CreateHorizontal_ConcatImage(imgs_cell)


% Concatenate images with spacers in between
img_sizes = cell2mat(cellfun(@(x) size(x),imgs_cell,'UniformOutput',0)');
bw_imgs = false(max(img_sizes,[],1) .* [1 2*numel(imgs_cell)-1]);
bw_rois = bw_imgs;

col = 1;
for n = 1:numel(imgs_cell)
    bw_imgs(1:img_sizes(n,1),col:col+img_sizes(n,2)-1) = imgs_cell{n};
     bw_rois(1:img_sizes(n,1),col:col+img_sizes(n,2)-1) = true(img_sizes(n,:)); 
    col = col+2*max(img_sizes(:,2));
end


end
