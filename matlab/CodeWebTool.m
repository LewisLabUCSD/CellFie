function[score, score_binary ,taskInfos]=CodeWebTool(data,ThreshType,SampleNumber,ref,percentile_or_value,percentile,value,EnoughSamples,LocalThresholdType,percentile_low,percentile_high,value_low,value_high)
% data.gene % entrez_ID associated to the gene
% data.value  % expression associated to the gene, genes x samples

% %Required 
% ThreshType: 'global' 
% SampleNumber: 7 
% 
% %Global Parameters 
% percentile_or_value: 'percentile'/'value'
% percentile: 50 
% value: 5 
% 
% %Local Parameters 
% EnoughSamples: 'no' 
% LocalThresholdType: 'MinMaxMean'/'Mean' 
% percentile_or_value: 'percentile'/'value'
% percentile_low: 50 
% percentile_high: 50 
% value_low: 5 
% value_high: 5 
% %ref: 'MT_recon_2_2_entrez.mat'

%load the info about the task structure
load('taskStructure')
taskInfos=struct2cell(taskStructure);
taskInfos=taskInfos';
taskInfos(:,5:end)=[];


%load the reference model in matlab
load(strcat('resources/libraries/',ref));

%% Depending on the model, load the list of reactions associated with the task
% All these files will al the following format
% essentialRxnsbyTask_name_of_model
%load(strcat('essentialRxnsbyTask_recon2_2')
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

% introduce a warning in the webtool about how many genes provided are
% mapped to the model
if isempty(gene_notInModel)
display('All genes provided are included in the reference model')
else
    display([num2str(length(gene_notInModel)),' genes provided are not included in the reference model:'])
    data.gene(gene_notInModel)
end

if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end

% remove the gene and associated data provided by the user that are not in the model
data.gene(gene_notInModel)=[];
if SampleNumber==1
    data.value(gene_notInModel)=[];
else
    data.value(gene_notInModel,:)=[];
end

if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end

% If there are negative values, shift positive
if any(any(data.value <= 0))
    warning('there are negative expression values, adding mininmum + 1e-2 to all values')
    data.value = data.value + abs(min(min(data.value))) + 1e-2;
end

% get the threshold value and the histogram for the complete dataset and
% print a figure
if SampleNumber>1
    linData = reshape(data.value,numel(data.value),1);
else
    linData=data.value;
end
linData(linData==0)=[];
LocalThresholdType
percentile_or_value
ThreshType
% definition of the thresholds
if strcmp(ThreshType,'global') && strcmp(percentile_or_value,'percentile') 
    display('RUN - global: percentile')
    figure(1);hist(log10(linData),50)
    l_global = (prctile(log10(linData),percentile));
    data.ths=10^l_global;
    h_global=line([l_global,l_global],get(gca,'YLim'),'Color','c');
    legend(h_global,{[num2str(percentile),'th percentile = ',num2str(data.ths)]},'FontSize',14)
    xlabel('log10(expressionValue)','FontSize',14)
    ylabel('Genes','FontSize',14)
elseif strcmp(ThreshType,'global') && strcmp(percentile_or_value,'value') 
    display('RUN - global: value')
    figure(1);hist(log10(linData),50)
    l_global = (prctile(log10(linData),percentile));
    data.ths=value;
    h_global=line([log10(value),log10(value)],get(gca,'YLim'),'Color','c');
    legend(h_global,{'global threshold'},'FontSize',14)
    xlabel('log10(expressionValue)','FontSize',14)
    ylabel('Genes','FontSize',14)
elseif strcmp(ThreshType,'local') && strcmp(LocalThresholdType,'Mean') 
    display('RUN - local: mean')
    figure(1);hist(log10(linData),50)
    xlabel('log10(expressionValue)','FontSize',14)
    ylabel('Genes','FontSize',14)
elseif strcmp(ThreshType,'local') && strcmp(LocalThresholdType,'minmaxmean')&& strcmp(percentile_or_value,'percentile') 
    display('RUN - local: minmaxmean: percentile')
    figure(1);hist(log10(linData),50)
    l_high = (prctile(log10(linData),percentile_high));
    data.ths_high=10^l_high;
    l_low = (prctile(log10(linData),percentile_low));
    data.ths_low=10^l_low;
    h_high=line([l_high,l_high],get(gca,'YLim'),'Color','c');
    h_low=line([l_low,l_low],get(gca,'YLim'),'Color','r');
    legend([h_high,h_low],{[num2str(percentile_high),'th percentile = ',num2str(data.ths_high)],[num2str(percentile_low),'th percentile = ',num2str(data.ths_low)]},'FontSize',14)
    xlabel('log10(expressionValue)','FontSize',14)
    ylabel('Genes','FontSize',14)
elseif strcmp(ThreshType,'local') && strcmp(LocalThresholdType,'minmaxmean')&& strcmp(percentile_or_value,'value')
    display('RUN - local: minmaxmean: value') 
    figure(1);hist(log10(linData),50)
    data.ths_high=value_high;
    data.ths_low=value_low;
    h_high=line([log10(value_high),log10(value_high)],get(gca,'YLim'),'Color','c');
    h_low=line([log10(value_low),log10(value_low)],get(gca,'YLim'),'Color','r');
    legend([h_high,h_low],{'Lower local threshold','Upper local threshold'},'FontSize',14)
    xlabel('log10(expressionValue)','FontSize',14)
    ylabel('Genes','FontSize',14)
else
    error('No analysis triggered')
end

saveas(figure(1),'Analysis/Figures/histogram.png')

if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end

%% Compute the threshold(s) depending on the approach used
Gene_score=[];
switch ThreshType
    case 'local'
        if strcmp(LocalThresholdType,'Mean')
            threshold=[];
            for i=1:length(data.gene) 
            	threshold(i)=mean(expressionValue);
            end
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

if size(data.value,1) ~= length(data.gene)
    error('data.value does not have the same number of rows as data.gene')
end

% Mapping of the expression data to the model
expression.gene=data.gene;
expression.Rxns=[];
expression.gene_used=[];
expression.count=[];
parsedGPR = GPRparser(model);% Extracting GPR data from model
for i=1:SampleNumber
    if SampleNumber==1
        expression.value=Gene_score;
    else
       expression.value=Gene_score(:,i);
    end
    [expressionRxns, gene_used] = mapExpressionToReactions_outParse(model,expression,parsedGPR); 
    gene_all=[];
    for j=1:length(gene_used)
        if ~isempty(gene_used{j})
           gene_all{end+1}=gene_used{j};
        end
    end
 %   unique_gene=unique(gene_all);
 %   unique_gene=unique_gene;
    %countGene = histc(gene_all(:), unique_gene);
    countGene = tabulate(gene_all);
    count=[];
    for k=1:length(gene_used)
        if ~isempty(gene_used{k})
           %count(k)=countGene(find(unique_gene==gene_used{k})); % removed
           %and replaced with line 190
           tmp=countGene(strcmp(countGene(:,1),gene_used{k}),2);
           count(k)=tmp{1};
        else
          count(k)=0;
        end
    end
    expression.Rxns=[expression.Rxns expressionRxns];
    expression.gene_used=[expression.gene_used gene_used];
    expression.count=[expression.count count'];
end

%% Compute the score
expressionRxns=expression.Rxns;
significance=1./expression.count;
significance(isinf(significance))=0;
ScorebyTask=[];
ScorebyTask_binary=[];
for i=1:210
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
noTask=find(ScorebyTask(:,1)==-1);
score=ScorebyTask;
score_binary=ScorebyTask_binary;
score([noTask],:) = [];
score_binary([noTask],:) = [];
taskInfos([noTask],:)=[];


end

