# RUN: 
# cd genepattern/
# ./bashwrapper.sh 'test/suite/dataTest.mat' 3 'test/suite/MT_recon_2_2_entrez.mat' 'local' 'value' 'minmaxmean' 25 75
## read parameters
#load('test/suite/dataTest.mat')
DATA=$1
#SampleNumber=3;
SAMP=$2
#ref='test/suite/MT_recon_2_2_entrez.mat';
REF=$3
#param.ThreshType='local';
pTHRESH=$4
#param.percentile_or_value='value';
pPERCVAL=$5
#param.LocalThresholdType='minmaxmean';
pTYPE=$6
#param.value_low=25;
pLOW=$7
#param.value_high=75;
pHIGH=$8

#Run
CMD="addpath(genpath('..'));execCellfie('$DATA',$SAMP,'$REF','$pTHRESH','$pPERCVAL','$pTYPE',$pLOW,$pHIGH);exit"
echo $CMD
matlab -nosplash -nodesktop -r $CMD


