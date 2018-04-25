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
[score, score_binary ,taskInfos]=CodeWebTool(data,ThreshType,SampleNumber,ref,percentile_or_value,percentile,value,EnoughSamples,LocalThresholdType,percentile_low,percentile_high,value_low,value_high)

rownames=taskInfos(:,2)
colnames=1:SampleNumber
score_long=wide2long(score,rownames,colnames)


save 'Analysis/Output/score.mat' score
%xlswrite('Analysis/Output/score.xlsx',score)
%csvwrite('Analysis/Output/score.csv',score)
%csvwrite('Analysis/Output/score.long.csv',score_long,{'Sample','Task','Value'})

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
% ax.YTickLabelRotation=90;

ax.XTick=1:length(colnames);
ax.XTickLabel=colnames(order3);
colorbar
saveas(figure(2),'Analysis/Figures/score_heatmap.png')

figure(3)
dendrogram(Z3)
saveas(figure(3),'Analysis/Figures/score_dendro.png')

%figure(4)
%CGobj = clustergram(score,...
%    'RowLabels', rownames,'ColumnLabels', colnames); %, ...
    %'RowPDist', RowPDistValue, 'ColumnPDist', ColumnPDistValue)
%saveas(figure(4),'Analysis/Figures/score_clustergram.png')

save 'Analysis/Output/score_binary.mat' score_binary

D2 = pdist(score_binary);
Z2 = linkage(D2,'average');
order2 = optimalleaforder(Z2, D2);
D3 = pdist(score_binary');
Z3 = linkage(D3,'average');
order3 = optimalleaforder(Z3, D3);

score_bin_cluster = score_binary(order2,order3);

figure(2)
imagesc(score_bin_cluster)
ax=gca;
ax.YTick=1:length(taskInfos(:,2));
ax.YTickLabel=taskInfos(order2,2);
% ax.YTickLabelRotation=90;

ax.XTick=1:length(colnames);
ax.XTickLabel=colnames(order3);
colorbar
saveas(figure(2),'Analysis/Figures/score_bin_heatmap.png')

figure(3)
dendrogram(Z3)
saveas(figure(3),'Analysis/Figures/score_bin_dendro.png')

%figure(4)
%CGobj = clustergram(score_binary,...
%    'RowLabels', rownames,'ColumnLabels', colnames); %, ...
%    %'RowPDist', RowPDistValue, 'ColumnPDist', ColumnPDistValue)
%saveas(figure(4),'Analysis/Figures/score_bin_clustergram.png')


%xlswrite('Analysis/Output/score_binary.xlsx',score_binary)
%csvwrite('Analysis/Output/score_binary.csv',score_binary)
%csvwrite('Analysis/Output/score_binary.long.csv',wide2long(score_binary,rownames,colnames),{'Sample','Task','Value'})
save 'Analysis/Output/taskInfos.mat' taskInfos
%xlswrite('Analysis/Output/taskInfos.xlsx',taskInfos)
%csvwrite('Analysis/Output/taskInfos.csv',taskInfos)
