function []=execCellfie(DATA,SAMP,REF,pTHRESH,pPERCVAL,pTYPE,pLOW,pHIGH,outputdir)
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
    if strcmp(pTHRESH,'local')
        if strcmp(pPERCVAL,'percentile')
            param.percentile_low=25;
            param.percentile_high=75;
        elseif strcmp(pPERCVAL,'value')
            param.value_low=str2num(pLOW);
            param.value_high=str2num(pHIGH);
        else
            error("cutoff type must be 'percentile' or 'value'")
        end
    elseif strcmp(pTHRESH,'global')
        if strcmp(pPERCVAL,'percentile')
            param.percentile=str2num(pLOW);
        elseif strcmp(pPERCVAL,'value')
            param.value=str2num(pLOW);
        else
            error("cutoff type must be 'percentile' or 'value'")
        end
    else
        error("threshold type must be 'local' or 'global'")
    end
	[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param);

	save cellfieout score score_binary taskInfos detailScoring
	%saveas(figure(1),'histogram.png')
	%close(figure(1))

	csvwrite(strcat(outputdir,'/score.csv'),score);
	csvwrite(strcat(outputdir,'/score_binary.csv'),score_binary);
	T = cell2table(taskInfos);
 	writetable(T,strcat(outputdir,'/taskInfo.csv'));

% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh \ 
%   /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.mat 3 \
%   MT_recon_2_2_entrez.mat local value minmaxmean 25 75 outtmp
% execCellfie('test/suite/dataTest.mat','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75','outtmp')
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh \
%   /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.csv 3 \
%   MT_recon_2_2_entrez.mat local value minmaxmean 25 75 outtmp
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75','outtmp')
