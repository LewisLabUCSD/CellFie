function []=runCellfie(DATA,SAMP,REF,pTHRESH,pPERCVAL,pTYPE,pLOW,pHIGH)
	%% dataRecon22_local_minmaxmean_value
	load(DATA)
	SampleNumber=SAMP;
	ref=REF;
	param.ThreshType=pTHRESH;
	param.percentile_or_value=pPERCVAL;
	param.LocalThresholdType=pTYPE;
	param.value_low=pLOW;
	param.value_high=pHIGH;

	[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);

	csvwrite('score.csv',score);
	csvwrite('score_binary.csv',score_binary);
	csvwrite('taskInfo.csv',taskInfos);
