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
