# CellFie

## Install

- Install matlab 2016b

- Install cobra

Download git,curl and cobratoolbox
```bash
# bash
sudo apt-get install git-all curl
git clone https://github.com/opencobra/cobratoolbox.git <desired path to cobra>/cobratoolbox
```
run initialization in matlab
```matlab
% matlab
cd <path to cobra>
initCobraToolbox
addpath(genpath(<path to cobra>));
savepath
```
Note: installing the solvers is not necessary for cellfie

- Install cellfie
Download Cellfie
```bash
# bash
git clone https://github.com/ResearchSoftwareInstitute/CellFie.git <desired path to cellfie>/CellFie
```
Initialize Cellfie
```matlab
% matlab
cd <path to cellfie>
initCellFie
addpath(genpath(<path to cellfie>));
savepath
```
## Quick Start
Run in matlab
```matlab
#matlab
load('test/suite/dataTest.mat')
SampleNumber=3;
ref='test/suite/MT_recon_2_2_entrez.mat';
param.ThreshType='local';
param.percentile_or_value='value';
param.LocalThresholdType='minmaxmean';
param.value_low=25;
param.value_high=75;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
```
Run with bash wrapper
```bash
#bash
cd genepattern/
./bashwrapper.sh 'dataTest.mat' 3 'MT_recon_2_2_entrez.mat' 'local' 'value' 'minmaxmean' 25 75
```
