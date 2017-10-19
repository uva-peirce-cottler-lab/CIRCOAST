function labelCCMatrix(cc, num_label)

bw_seg = false(cc.ImageSize);
for n=(1:cc.NumObjects)
    bw_seg(cc.PixelIdxList{n})=1;
end
rpstat_st = regionprops(bw_seg,'Centroid');
%         hole_area = arrayfun(@(x) x.Area(1),rpstat_st);
centroid = [arrayfun(@(x) x.Centroid(1),rpstat_st) arrayfun(@(x) x.Centroid(2),rpstat_st)];


%         [A,ix] = sort(hole_area);
%         figure; imshow(~bw_skel);
figure
imshow(bw_seg) 
hold on
text(centroid(:,1),centroid(:,2),...
    arrayfun(@(x) sprintf('%i', x),num_label,'UniformOutput',0),...
    'HorizontalAlignment','center','VerticalAlignment','middle','Color','red');
%         text(centroid(ix,1),centroid(ix,2),arrayfun(@(x) sprintf('%i',x),hole_area(ix),'UniformOutput',0),...
%             'HorizontalAlignment','center','VerticalAlignment','bottom','Color','black')
%
hold off

