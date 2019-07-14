# CellFie (FOR GENE PATTERN)

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
%matlab
% expression matrix: entrez ids x samples
load 'test/suite/dataTest.mat'
% number of samples (equal to the column number of the expression matrix
SampleNumber=3;

% reference genome (all listed in the test/suite)
ref='test/suite/MT_recon_2_2_entrez.mat';

% type of thresholding method ('local' or 'global')
param.ThreshType='local';
% local threshold type ('minmaxmean', 'custom' or 'mean')
param.LocalThresholdType='minmaxmean';

% numerical cutoff ('value' or 'percent')
param.percentile_or_value='value';
% minimum cutoff (0-1000 for 'value', 0-1 'percent')
param.value_low=25;
% maximum cutoff (0-1000 for 'value', 0-1 'percent')
param.value_high=75;

[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);
```
Run with bash wrapper
```bash
#bash
cd genepattern/
./bashwrapper.sh 'dataTest.mat' 3 'MT_recon_2_2_entrez.mat' 'local' 'value' 'minmaxmean' 25 75
```
## [Explanation of method and parameters (wiki)](https://github.com/ResearchSoftwareInstitute/CellFie/wiki/Cellfie-Documentation)
