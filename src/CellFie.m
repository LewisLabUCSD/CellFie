function[score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param)
% Compute the score associated to each metabolic task listed in taskstructure based on transcriptomic data
%
% USAGE:
%    [score, score_binary ,taskInfos, detailScoring]=CellFie(data,SampleNumber,ref,param)
%
% INPUTS:
%	data
%       .gene                   cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  mRNA expression data structure (genes x samples)associated to each gene metioned in data.gene
%   SampleNumber                Number of samples
%   ref                         Reference model used to compute the
%                               metabolic task scores (e.g.,'MT_recon_2_2_entrez.mat')
% OPTIONAL INPUTS:
%   param.ThreshType            Type of thresholding approach used
%                               (i.e.,'global' or 'local') (default - local)
% related to the use of a GLOBAL thresholding approach - the threshold value is the same for all the genes
%   param.percentile_or_value   the threshold can be defined using a value introduced by the user ('value')
%                               or based on a percentile of the distribution of expression value for all the
%                               genes and across all samples of your
%                               dataset ('percentile')
%   param.percentile            percentile from the distribution of
%                               expression values for all the genes and across all samples that will be
%                               used to define the threshold value
%   param.value                 expression value for which a gene is
%                               considered as active or not (e.g., 5)
%
% related to the use of a LOCAL thresholding approach - the threshold value is different for all the genes
%   param.percentile_or_value   the threshold can be defined using a value introduced by the user ('value')
%                               or based on a percentile of the distribution of expression value of a 
%                               specific gene across all samples of your
%                               dataset ('percentile'-default)
%   param.LocalThresholdType    option to define the type of local thresholding approach to use
%                               - 'minmaxmean' (default options )- the threshold for a gene is determined by the mean of expression
%                                  values observed for that gene among all the samples, tissues, or conditions BUT
%                                   the threshold :(i) must be higher or equal to a lower bound and (ii) must be lower
%                                   or equal to an upper bound.
%                               - 'mean' -the threshold for a gene is defined as the mean expression value
%                                  of this gene across all the samples, tissues, or conditions
%   param.percentile_low        lower percentile used to define which gene
%                               are always inactive in the case of use 'MinMaxMean' local thresholding
%                               approach (default = 25)
%   param.percentile_high       upper percentile used to define which gene
%                               are always active in the case of use 'MinMaxMean' local thresholding
%                               approach (default= 75)
%   param.value_low             lower expression value used to define which gene
%                               are always inactive in the case of use 'MinMaxMean' local thresholding
%                               approach (e.g., 5)
%   param.value_high            upper expression value used to define which gene
%                               are always active in the case of use 'MinMaxMean' local thresholding
%                               approach (e.g., 5)
%
% related to the gene mapping approach used
%   param.minSum:               instead of using min and max, use min for AND and Sum
%                               for OR (default: false, i.e. use min)
% OUTPUTS:
%   score                       relative quantification of the activity of a metabolic task in a specific condition
%                               based on the availability of data for multiple conditions
%   score_binary                binary version of the metabolic task score
%                               to determine whether a task is active or inactive in specific
%                               conditions
%   taskInfos                   Description of the metabolic task assessed
%   detailScoring               Matrix detailing the scoring
%       1st column = sample ID
%       2nd column = task ID
%       3th column = task score for this sample
%       4th column = task score in binary version for this sample
%       5th column = essential reaction associated to this task
%       6th column = expression score associated  to the reaction listed in the 5th column
%       7th column = gene used to determine the expression of the reaction listed in the 5th column
%       8th column = original expression value of the gene listed in the 7th column
%
% .. Authors:
%       - Anne Richelle, January 2019
if size(data.value,2) ~= SampleNumber
    error('The number of samples defined is not the same as the size of the dataset')
end
if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end
if ~exist('ref','var')
    error('The reference model has not been defined - please choose a reference model')
end
if ~exist('param','var')
    param.ThreshType='local';
    param.percentile_or_value='percentile';
    param.LocalThresholdType='minmaxmean';
    param.percentile_low=25;
    param.percentile_high=75;
%     percentile=[];
%     value=[];
%     value_low=[];
%     value_high=[];
%     minSum = false;
end

%load the info about the task structure
load('taskStructure')
taskInfos=struct2cell(taskStructure);
taskInfos=taskInfos';
taskInfos(:,5:end)=[];


%load the reference model in matlab
%%%% TO DO !!! DEFINE the good path to search the model
load(ref);

%% Depending on the model, load the list of reactions associated with the task
% All these files have the following format
% essentialRxnsbyTask_name_of_model and are located in essentialRxns folder
load(strcat('essentialRxns/essentialRxnsbyTask_',ref));

% check that at least part of list of gene provided are in the model loaded
ID_model=[];
gene_notInModel=[];
for i=1:length(data.gene)
    if isempty(find(strcmp(data.gene{i},model.genes)))
        gene_notInModel(end+1)=i;
    else
        tmpid=find(strcmp(data.gene{i},model.genes)==1);
        ID_model(end+1)=tmpid(1);
    end
end

% introduce a warning in the webtool about how many of the genes provided are
% actually mapped to the model
if isempty(gene_notInModel)
    display('All genes provided in data are included in the reference model')
else
    display([num2str(length(gene_notInModel)),' genes provided are not included in the reference model:'])
    data.gene(gene_notInModel)
end

% remove the gene and associated value provided by the user in data that are not in the model
data.gene(gene_notInModel)=[];
if SampleNumber==1
    data.value(gene_notInModel)=[];
else
    data.value(gene_notInModel,:)=[];
end

% get the threshold value and the histogram for the complete dataset and
% print a figure
if SampleNumber>1
    linData = reshape(data.value,numel(data.value),1);
else
    linData=data.value;
end
linData(linData==0)=[];

% definition of the thresholds
if strcmp(param.ThreshType,'global') && strcmp(param.percentile_or_value,'percentile')
    display('RUN - global: percentile')
    %figure(1);hist(log10(linData),50)
    l_global = (prctile(log10(linData),param.percentile));
    data.ths=10^l_global;
    %h_global=line([l_global,l_global],get(gca,'YLim'),'Color','c');
    %legend(h_global,{[num2str(param.percentile),'th percentile = ',num2str(data.ths)]},'FontSize',14)
    %xlabel('log10(expressionValue)','FontSize',14)
    %ylabel('Genes','FontSize',14)
elseif strcmp(param.ThreshType,'global') && strcmp(param.percentile_or_value,'value')
    display('RUN - global: value')
    %figure(1);hist(log10(linData),50)
    data.ths=param.value;
    %h_global=line([log10(param.value),log10(param.value)],get(gca,'YLim'),'Color','c');
    %legend(h_global,{'global threshold'},'FontSize',14)
    %xlabel('log10(expressionValue)','FontSize',14)
    %ylabel('Genes','FontSize',14)
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'mean')
    display('RUN - local: mean')
    %figure(1);hist(log10(linData),50)
    %xlabel('log10(expressionValue)','FontSize',14)
    %ylabel('Genes','FontSize',14)
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'minmaxmean')&& strcmp(param.percentile_or_value,'percentile')
    display('RUN - local: minmaxmean: percentile')
    %figure(1);hist(log10(linData),50)
    l_high = (prctile(log10(linData),param.percentile_high));
    data.ths_high=10^l_high;
    l_low = (prctile(log10(linData),param.percentile_low));
    data.ths_low=10^l_low;
    %h_high=line([l_high,l_high],get(gca,'YLim'),'Color','c');
    %h_low=line([l_low,l_low],get(gca,'YLim'),'Color','r');
    %legend([h_high,h_low],{[num2str(param.percentile_high),'th percentile = ',num2str(data.ths_high)],[num2str(param.percentile_low),'th percentile = ',num2str(data.ths_low)]},'FontSize',14)
    %xlabel('log10(expressionValue)','FontSize',14)
    %ylabel('Genes','FontSize',14)
elseif strcmp(param.ThreshType,'local') && strcmp(param.LocalThresholdType,'minmaxmean')&& strcmp(param.percentile_or_value,'value')
    display('RUN - local: minmaxmean: value')
    %figure(1);hist(log10(linData),50)
    data.ths_high=param.value_high;
    data.ths_low=param.value_low;
    %h_high=line([log10(param.value_high),log10(param.value_high)],get(gca,'YLim'),'Color','c');
    %h_low=line([log10(param.value_low),log10(param.value_low)],get(gca,'YLim'),'Color','r');
    %legend([h_high,h_low],{'Lower local threshold','Upper local threshold'},'FontSize',14)
    %xlabel('log10(expressionValue)','FontSize',14)
    %ylabel('Genes','FontSize',14)
else
    error('No analysis triggered')
end

%% TO DO - DEFINE THE GOOD PATH TO SAVE FIGURE -
% saveas(figure(1),'Analysis/Figures/histogram.png')

%% Compute the threshold(s) depending on the approach used
Gene_score=[];
switch param.ThreshType
    case 'local'
        if strcmp(param.LocalThresholdType,'mean')
            %the threshold for each gene is equal to its mean value over
            %all the samples
            threshold=mean(data.value,2)';
        else
            threshold=[];
            for i=1:length(data.gene)
            	expressionValue=data.value(i,:);
               	if  mean(expressionValue)>=data.ths_high
                	threshold(i)=data.ths_high;
                else
                	threshold(i)=max(mean(expressionValue),data.ths_low);
             	end
            end
        end
    % every single gene is associated to an expression score
    for i=1:SampleNumber
        Gene_score(:,i)=5.*log(1+(data.value(:,i)./threshold'));
    end
    case 'global'
        Gene_score=5.*log(1+(data.value./data.ths));
end

% Mapping of the expression data to the model
expression.gene=data.gene;
expression.Rxns=[];
expression.gene_used=[];
expression.count=[];
minSum = false;

display('Load GPR parse');
parsedGPR = GPRparser(model,minSum);% Extracting GPR data from model % Ben Kellman, 8/25/19

display('Mapping of the expression data to the model');
for i=1:SampleNumber
    if SampleNumber==1
        expression.value=Gene_score;
    else
       expression.value=Gene_score(:,i);
    end
    %[expressionRxns, parsedGPR, gene_used] =
    %mapExpressionToReactions(model, expression, minSum); % moved gpr
    %parsing outside of for loop $ Ben Kellman, 8/25/19
    % moved out of forloop
    % Find wich genes in expression data are used in the model
    [gene_id, gene_expr] = findUsedGenesLevels(model,expression);
    % Link the gene to the model reactions # NOTE: THIS IS THE TIME SINK
    [expressionRxns, gene_used] = selectGeneFromGPR(model, gene_id, gene_expr, parsedGPR, minSum);
    
    gene_all=[];
    for j=1:length(gene_used)
        if ~isempty(gene_used{j})
           gene_all(end+1)=str2num(gene_used{j}{1});
        end
    end
    countGene = tabulate(gene_all);
    count=[];
    for k=1:length(gene_used)
        if ~isempty(gene_used{k})
           tmp=countGene(str2num(gene_used{k}{1}),2);
           count(k)=tmp;
        else
          count(k)=0;
        end
    end
    expression.Rxns=[expression.Rxns expressionRxns];
    expression.gene_used=[expression.gene_used gene_used'];
    expression.count=[expression.count count'];
end

%% Compute the score
expressionRxns=expression.Rxns;
significance=1./expression.count;
significance(isinf(significance))=0;
ScorebyTask=[];
ScorebyTask_binary=[];
display('Compute the task activity score');
for i=1:size(taskInfos,1)
	if ~isempty(essentialRxns{i})
    	rxns=essentialRxns{i};
        rxnID=findRxnIDs(model,rxns);
    	rxnID(rxnID==0)=[];
      	if ~isempty(rxnID)
        	expValue=expressionRxns(rxnID,:);
           	signValue=significance(rxnID,:);
            % if no gene is associated with one of the reaction -
          	% remove the reactions from the count
            if ~isempty(find(sum(expValue,2)==-SampleNumber))
                signValue(find(sum(expValue,2)==-SampleNumber),:)=[];
                expValue(find(sum(expValue,2)==-SampleNumber),:)=[];
            end
           	if ~isempty(expValue)
            	if size(expValue,1)>1
                    ScorebyTask(i,:)=sum(expValue.*signValue)./size(expValue,1);
                    Val=sum(expValue)./size(expValue,1);
                    ID_up=find(Val>=5*log(2));
                    ScorebyTask_binary(i,:)=zeros(1,SampleNumber);
                    ScorebyTask_binary(i,ID_up)=1;
                else
                    ScorebyTask(i,:)=expValue.*signValue;
                    ID_up=find(expValue>=5*log(2));
                    ScorebyTask_binary(i,:)=zeros(1,SampleNumber);
                    ScorebyTask_binary(i,ID_up)=1;
                end
            else
            	ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
                ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
           	end
        else
        	ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
            ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
    	end
    else
        ScorebyTask(i,:)=-1.*ones(1,SampleNumber);
        ScorebyTask_binary(i,:)=-1.*ones(1,SampleNumber);
   end
end

detailScoring={};
display('Format the score')
for j=1:SampleNumber
    incR=1;
    for i=1:size(taskInfos,1)
        if ~isempty(essentialRxns{i})
            rxns=essentialRxns{i};
            rxnID=findRxnIDs(model,rxns);
            rxnID(rxnID==0)=[];
            if ~isempty(rxnID)
                for k=1:length(rxnID)
                    %1st column = sample ID
                    detailScoring{incR,((j-1)*8)+1}=j;
                    %2nd column = task ID
                    detailScoring{incR,((j-1)*8)+2}=i;
                    %3th column = task score for this sample
                    detailScoring{incR,((j-1)*8)+3}=ScorebyTask(i,j);
                    %4th column = task score in binary version for this sample
                    detailScoring{incR,((j-1)*8)+4}=ScorebyTask_binary(i,j);
                    %5th column = essential reaction associated to this
                    %task
                    detailScoring{incR,((j-1)*8)+5}=rxns(k);
                    %6th column = expression score associated  to the
                    %reaction listed in the 5th column
                    detailScoring{incR,((j-1)*8)+6}=expression.Rxns(rxnID(k),SampleNumber);
                    %7th column = gene used to determine the expression of the
                    %reaction listed in the 5th column
                    geneName=expression.gene_used(rxnID(k),SampleNumber);
                    detailScoring{incR,((j-1)*8)+7}=geneName;
                    %8th column = original expression value of the gene
                    %listed in the 7th column
                    detailScoring{incR,((j-1)*8)+8}=data.value(strcmp(data.gene,geneName),SampleNumber);
                    incR=incR+1;
                end
           end
        end
    end
    %noTask=find(ScorebyTask(:,1)==-1);
    score=ScorebyTask;
    score_binary=ScorebyTask_binary;
    %score([noTask],:) = [];
    %score_binary([noTask],:) = [];
    %taskInfos([noTask],:) = [];
end

% get model checking tasks
r = readtable('input/inactive_reactions.csv','Delimiter',',');
noTask=r.index; %find(ScorebyTask(:,1)==-1);
% check removed tasks
for i=1:length(noTask)
    nti = noTask(i);
    if taskInfos{nti,2}~=r.v1{i}
        error('inactive reaction index does not match taskInfo')
    end
end
    
% remove model checking tasks
score([noTask],:) = [];
score_binary([noTask],:) = [];
taskInfos([noTask],:) = [];
