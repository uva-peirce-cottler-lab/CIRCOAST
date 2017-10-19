clear all

ntrials=500;
r=100;
hdat = zeros(r,10);
bino_cdf = binocdf(0:ntrials,ntrials,.5);

for n=1:r
    
    
    bino_numbers = binornd(ntrials,.5,[1 1000000]);
    
    p_vals = 1-bino_cdf(bino_numbers(1:100)+1);
    
    hdat(n,:) = hist(p_vals,10);
end

save data.mat

load data.mat

h = histogram(p_vals,10);
% h.Values = mean(hdat,1);

bar_edges = h.BinEdges;

hb = bar(bar_edges, [mean(hdat,1) 0]./sum(mean(hdat,1)),'histc')
set(hb, 'FaceColor', [102 170 215]/255)
xlabel('p Value')
ylabel('Probabbility')
axis([0 1 0 .15])

