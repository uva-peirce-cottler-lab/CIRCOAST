function cell_overlap_chooser()


bw = false(300,600);

cell_rad = 100;
cell_rat = [3:-.1:0];

for n = 1:numel(cell_rat)
% Add first cell
bw(cell_rad*1.5,cell_rad*1.5)=1;
bw(cell_rad*1.5,cell_rat(n)*cell_rad+cell_rad*1.5)=1;
bw = bwdist(bw)<=100;

% bw = imdilate(bw,strel('disk', cell_rad));
% imshow(bw)

imwrite(bw,['radius_' sprintf('%2.1f',cell_rat(n)) '.tif']);
bw(:)=0;
end