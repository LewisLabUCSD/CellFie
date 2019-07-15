function []=runCellfie(DATA,SAMP,REF,pTHRESH,pPERCVAL,pTYPE,pLOW,pHIGH)
	load(DATA);
	SampleNumber=str2num(SAMP);
	ref=REF;
	param.ThreshType=pTHRESH;
	param.percentile_or_value=pPERCVAL;
	param.LocalThresholdType=pTYPE;
	param.value_low=str2num(pLOW);
	param.value_high=str2num(pHIGH);

	[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);

	%save cellfieout score score_binary taskInfos detailScoring
	saveas(figure(1),'histogram.png')
	close(figure(1))

	csvwrite('score.csv',score);
	csvwrite('score_binary.csv',score_binary);
	T = cell2table(taskInfos);
 	writetable(T,'taskInfo.csv');

% ./genepattern/CellFie/for_redistribution_files_only/run_runCellFie.sh /media/ben/9c17f1c9-a45e-49ec-b547-8fbd2f25ccc6/tmp/v94 test/suite/dataTest.mat 3 MT_recon_2_2_entrez.mat local value minmaxmean 25 75
