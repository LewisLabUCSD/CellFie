addpath(genpath('/var/www/pinapl-py-site/'))
addpath(genpath('/var/www/MATLAB/lib/'))
addpath(genpath('/var/www/MATLAB/lib/MetTasks/'))
addpath(genpath('/var/www/MATLAB/lib/MetTasks/matlab'))
addpath(genpath('/var/www/MATLAB/lib//MetTasks/matlab/yaml'))
addpath(genpath('/var/www/MATLAB/lib/cobratoolbox/'))


yml = YAML.read(config_file);

dataAll = readtable(data_file);
'input data head:'
dataAll(1:10,:)

ThreshType = yml.ThreshType
SampleNumber = yml.SampleNumber
ref = yml.ref

percentile_or_value = yml.percentile_or_value
percentile = yml.percentile
value = yml.value

EnoughSamples = yml.EnoughSamples
LocalThresholdType = yml.LocalThresholdType
percentile_low = yml.percentileLow
percentile_high = yml.percentileHigh
value_low = yml.valueLow
value_high = yml.valueHigh

data = {};
data.gene = string( table2array( dataAll(:,1) ) );
data.value = table2array( dataAll(:,2:end) );
[score, score_binary ,taskInfos, detailScoring]=CellFie(data,ThreshType,SampleNumber,ref,percentile_or_value,percentile,value,LocalThresholdType,percentile_low,percentile_high,value_low,value_high)

rownames=taskInfos(:,2)
colnames=1:SampleNumber
score_long=wide2long(score,rownames,colnames)

fname='Analysis/Output/'
save 'Analysis/Output/score.mat' score
tab=table(taskInfos,score)
writetable(tab,fullfile(fname, 'score.csv'))

% write names if present
tableNames = dataAll.Properties.VariableNames
writetable(tab,fullfile(fname, 'annotation.csv'))

D2 = pdist(score);
Z2 = linkage(D2,'average');
order2 = optimalleaforder(Z2, D2);
D3 = pdist(score');
Z3 = linkage(D3,'average');
order3 = optimalleaforder(Z3, D3);

score_cluster = score(order2,order3);

figure(2)
imagesc(score_cluster)
ax=gca;
ax.YTick=1:length(taskInfos(:,2));
ax.YTickLabel=taskInfos(order2,2);

ax.XTick=1:length(colnames);
ax.XTickLabel=colnames(order3);
colorbar
saveas(figure(2),'Analysis/Figures/score_heatmap.png')

figure(3)
dendrogram(Z3)
saveas(figure(3),'Analysis/Figures/score_dendro.png')

fname='Analysis/Output/'
save 'Analysis/Output/score_binary.mat' score_binary
tab=table(taskInfos,score_binary)
writetable(tab,fullfile(fname, 'score_binary.csv'))

D2 = pdist(score_binary);
Z2 = linkage(D2,'average');
order2 = optimalleaforder(Z2, D2);
D3 = pdist(score_binary');
Z3 = linkage(D3,'average');
order3 = optimalleaforder(Z3, D3);

score_bin_cluster = score_binary(order2,order3);

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(score_bin_cluster)
ax=gca;

ax.XTick=1:length(colnames);
ax.XTickLabel=colnames(order3);
colorbar
fig=gcf;
saveas(figure(fig.Number),'Analysis/Figures/score_bin_heatmap.png')

figure(3)
dendrogram(Z3)
saveas(figure(3),'Analysis/Figures/score_bin_dendro.png')

save 'Analysis/Output/taskInfos.mat' taskInfos
fname='Analysis/Output/'
tab=table(taskInfos)
writetable(tab,fullfile(fname, 'taskInfos.csv'))

exit
exit