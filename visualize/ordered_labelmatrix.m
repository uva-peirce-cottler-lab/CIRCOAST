function db = ordered_labelmatrix(cc)

db = zeros(cc.ImageSize);

for n=1:cc.NumObjects
   db(cc.PixelIdxList{n})=n; 
end
