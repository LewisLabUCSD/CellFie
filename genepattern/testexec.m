
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.mat 3 MT_recon_2_2_entrez.mat local value NA minmaxmean 25 75 outtmp
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.csv 3 MT_recon_2_2_entrez.mat local value NA minmaxmean 25 75 outtmp
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.xlsx 3 MT_recon_2_2_entrez.mat local value NA minmaxmean 25 75 outtmp
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.mat 3 MT_recon_2_2_entrez.mat local value NA minmaxmean 25 75 outtmp

execCellfie('../test/suite/dataTest.xlsx','3','MT_recon_2_2_entrez.mat','local','value','25','minmaxmean','25','75','outtmp')
execCellfie('../test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','percentile','40','minmaxmean','15','85','outtmp') 
execCellfie('../test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','value','40','mean','15','85','outtmp')
execCellfie('../test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','global','percentile','30','minmaxmean','15','85','outtmp')
