function setErrorBarWidth(hE)


for n=1:numel(hE)
% adjust error bar width
hE_c                   = ...
    get(hE(n), 'Children'    );

keyboard
errorbarXData          = ...
    get(hE_c, 'XData'       );
errorbarXData(4:9:end) = ...
    errorbarXData(1:9:end) - 0.2;
errorbarXData(7:9:end) = ....
    errorbarXData(1:9:end) - 0.2;
errorbarXData(5:9:end) = ...
    errorbarXData(1:9:end) + 0.2;
errorbarXData(8:9:end) = ...
    errorbarXData(1:9:end) + 0.2;
set(hE_c(2), 'XData', errorbarXData);


end