function zs = zscore_col(dat)

zs=zeros(size(dat));
for c=1:size(dat,2)
    zs(:,c)= zscore(dat(:,c));
end