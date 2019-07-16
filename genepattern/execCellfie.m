function []=execCellfie(DATA,SAMP,REF,pTHRESH,pPERCVAL,pTYPE,pLOW,pHIGH)
    if contains(DATA,'mat')
    	load(DATA);
    elseif( contains(DATA,'csv')||contains(DATA,'tsv'))
        datatmp=readtable(DATA);
        genetmp=table2array(datatmp(:,1));
        tmp = struct('gene',table2array(datatmp(:,1)),'value',table2array(datatmp(:,2:end)));
        tmp.gene = num2cell(tmp.gene);
        for i=1:length(genetmp)
            tmp.gene{i}=int2str(genetmp(i));
        end
        data=tmp;
    end
	SampleNumber=str2num(SAMP);
	ref=REF;
	param.ThreshType=pTHRESH;
	param.percentile_or_value=pPERCVAL;
	param.LocalThresholdType=pTYPE;
	param.value_low=str2num(pLOW);
	param.value_high=str2num(pHIGH);

	[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);

	save cellfieout score score_binary taskInfos detailScoring
	saveas(figure(1),'histogram.png')
	close(figure(1))

	csvwrite('score.csv',score);
	csvwrite('score_binary.csv',score_binary);
	T = cell2table(taskInfos);
 	writetable(T,'taskInfo.csv');

% ./genepattern/CellFie/for_redistribution_files_only/run_runCellFie.sh /media/ben/9c17f1c9-a45e-49ec-b547-8fbd2f25ccc6/tmp/v94 test/suite/dataTest.mat 3 MT_recon_2_2_entrez.mat local value minmaxmean 25 75
% execCellfie('test/suite/dataTest.mat','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75')
% ./genepattern/CellFie/for_redistribution_files_only/run_runCellFie.sh /media/ben/9c17f1c9-a45e-49ec-b547-8fbd2f25ccc6/tmp/v94 test/suite/dataTest.csv 3 MT_recon_2_2_entrez.mat local value minmaxmean 25 75
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75')
