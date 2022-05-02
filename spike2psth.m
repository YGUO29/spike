function psth = spike2psth(sp_train, window, nCluster, tend)
if isempty(tend)
    tend = max(sp_train + window);
end
edges = 0:window:tend;
[psth, edges] = histcounts(sp_train, edges);
psth = psth./window./nCluster;

figure('color','w','position',[105         659        2295         344])
% bar(edges(2:end), psth)
plot(edges(2:end), psth)
end

