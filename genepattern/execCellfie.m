function []=execCellfie(DATA,SAMP,REF,pTHRESH,pPERCVAL,pGLOBAL,pTYPE,pLOW,pHIGH,outputdir)
    if contains(DATA,'mat')
    	load(DATA);
    elseif( contains(DATA,'csv')||contains(DATA,'tsv')||contains(DATA,'xlsx')||contains(DATA,'xls'))
        datatmp=readtable(DATA);
        genetmp=table2array(datatmp(:,1));
        tmp = struct('gene',table2array(datatmp(:,1)),'value',table2array(datatmp(:,2:end)));
        tmp.gene = num2cell(tmp.gene);
        for i=1:length(genetmp)
            tmp.gene{i}=int2str(genetmp(i));
        end
        data=tmp;
    else
        error('Your data file must be formatted as: .mat, .csv, .tsv, .xlsx, or .xls with the correct suffix')
    end
	SampleNumber=str2num(SAMP);
	ref=REF;
	param.ThreshType=pTHRESH;
	param.percentile_or_value=pPERCVAL;
	param.LocalThresholdType=pTYPE;
    if strcmp(pTHRESH,'local')
    	if strcmp(pTYPE,'minmaxmean')
	    if strcmp(pPERCVAL,'percentile')
		param.percentile_low=str2num(pLOW);
		param.percentile_high=str2num(pHIGH);
	    elseif strcmp(pPERCVAL,'value')
		param.value_low=str2num(pLOW);
		param.value_high=str2num(pHIGH);
	    else
		error("cutoff type must be 'percentile' or 'value'")
	    end
	end
    elseif strcmp(pTHRESH,'global')
        if strcmp(pPERCVAL,'percentile')
            param.percentile=str2num(pGLOBAL);
        elseif strcmp(pPERCVAL,'value')
            param.value=str2num(pGLOBAL);
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
    D = cell2table(detailScoring);
 	writetable(T,strcat(outputdir,'/taskInfo.csv'));
 	writetable(D,strcat(outputdir,'/detailScoring.csv'));

% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh \ 
%   /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.mat 3 \
%   MT_recon_2_2_entrez.mat local value minmaxmean 25 75 outtmp
% execCellfie('test/suite/dataTest.mat','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75','outtmp')
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh \
%   /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.csv 3 \
%   MT_recon_2_2_entrez.mat local value minmaxmean 25 75 outtmp
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','value','minmaxmean','25','75','outtmp')
% ./matlab_compiled/execCellfie/for_redistribution_files_only/run_execCellfie.sh \
%   /usr/local/MATLAB/MATLAB_Runtime/v94 test/suite/dataTest.xlsx 3 \
%   MT_recon_2_2_entrez.mat local value minmaxmean 25 75 outtmp
% execCellfie('test/suite/dataTest.xlsx','3','MT_recon_2_2_entrez.mat','local','value','25','minmaxmean','25','75','outtmp')
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','percentile','40,'minmaxmean','15','85','outtmp')
% 
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','local','value','40','mean','15','85','outtmp')
% execCellfie('test/suite/dataTest.csv','3','MT_recon_2_2_entrez.mat','global','percentile','30','minmaxmean','15','85','outtmp')
