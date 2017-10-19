function pval = CELLCOAV_1S(x,PLOT)

elem = @(x,k) x(k);
num_permute = 1e7;

% Uniformly sample p values from 0-1 wotj same # samples as data
unif_p_sampling = unifrnd(0,1, [num_permute, numel(x)]);%-0.5;
unif_p_mean = mean(unif_p_sampling,2);
unif_p_pdf = sort(unif_p_mean);

% Figure out what percential actual p value is
bv = unif_p_pdf>(mean(x));
try
    ind = elem(elem(1:num_permute,bv),1);
catch ME;
    ind = num_permute; fprintf('ERROR:Exceed max index\n');
end
pval = ind/num_permute;

if PLOT
    figure;
    hist_data = histogram(unif_p_mean,50,'Normalization','probability');
    hold on
    y=ylim;
    plot([mean(x) mean(x)],...
        [0 1.2*y(2)],'r','LineWidth',2)
    xlabel(['Mean CELLCOAV p value'])
    ylabel('Probability')
    axis([0 1 1.2*y]);
    hold off
    
    set(findall(gcf,'-property','FontSize'),'FontSize',8);
    set(findall(gcf,'-property','FontName'),'FontName','Helvetica');
    set(gcf,'position', [100 100 260 175])
end

