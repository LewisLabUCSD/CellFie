% Generation of all the potential results

%% dataRecon22_local_minmaxmean_value
load('dataTest.mat')
SampleNumber=3;
ref='MT_recon_2_2_entrez.mat';
param.ThreshType='local';
param.percentile_or_value='value';
param.LocalThresholdType='minmaxmean';
param.value_low=25;
param.value_high=75;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
save dataRecon22_local_minmaxmean_value score score_binary taskInfos detailScoring
csvwrite('dataRecon22_local_minmaxmean_value.score.csv',score);
csvwrite('dataRecon22_local_minmaxmean_value.score_binary.csv',score_binary);
T = cell2table(taskInfos);
writetable(T,'dataRecon22_local_minmaxmean_value.taskInfo.csv');

%% dataRecon22_local_minmaxmean_percentile
load('dataTest.mat')
SampleNumber=3;
ref='MT_recon_2_2_entrez.mat';
param.ThreshType='local';
param.percentile_or_value='percentile';
param.LocalThresholdType='minmaxmean';
param.percentile_low=25;
param.percentile_high=75;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
save dataRecon22_local_minmaxmean_percentile score score_binary taskInfos detailScoring
csvwrite('dataRecon22_local_minmaxmean_percentile.score.csv',score);
csvwrite('dataRecon22_local_minmaxmean_percentile.score_binary.csv',score_binary);
T = cell2table(taskInfos);
writetable(T,'dataRecon22_local_minmaxmean_percentile.taskInfo.csv');


%% dataRecon22_local_mean
load('dataTest.mat')
SampleNumber=3;
ref='MT_recon_2_2_entrez.mat';
param.ThreshType='local';
param.LocalThresholdType='mean';

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
save dataRecon22_local_mean score score_binary taskInfos detailScoring
csvwrite('dataRecon22_local_mean.score.csv',score);
csvwrite('dataRecon22_local_mean.score_binary.csv',score_binary);
T = cell2table(taskInfos);
writetable(T,'dataRecon22_local_mean.taskInfo.csv');


%% dataRecon22_global_value
load('dataTest.mat')
SampleNumber=3;
ref='MT_recon_2_2_entrez.mat';
param.ThreshType='global';
param.percentile_or_value='value';
param.value=50;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
save dataRecon22_global_value score score_binary taskInfos detailScoring
csvwrite('dataRecon22_global_value.score.csv',score);
csvwrite('dataRecon22_global_value.score_binary.csv',score_binary);
T = cell2table(taskInfos);
writetable(T,'dataRecon22_global_value.taskInfo.csv');


%% dataRecon22_global_percentile
load('dataTest.mat')
SampleNumber=3;
ref='MT_recon_2_2_entrez.mat';
param.ThreshType='global';
param.percentile_or_value='percentile';
param.percentile=50;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
save dataRecon22_global_percentile score score_binary taskInfos detailScoring
csvwrite('dataRecon22_global_percentile.score.csv',score);
csvwrite('dataRecon22_global_percentile.score_binary.csv',score_binary);
T = cell2table(taskInfos);
writetable(T,'dataRecon22_global_percentile.taskInfo.csv');
