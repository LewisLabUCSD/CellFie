function [gene_id, gene_expr] = findUsedGenesLevels_modif(model, exprData, printLevel)
% Returns vectors of gene identifiers and corresponding gene expression
% levels for each gene present in the model ('model.genes').
%
% USAGE:
%    [gene_id, gene_expr] = findUsedGenesLevels(model, exprData)
%
% INPUTS:
%
%   model:               input model (COBRA model structure)
%
%   exprData:            mRNA expression data structure
%       .gene                cell array containing GeneIDs in the same
%                            format as model.genes
%       .value               Vector containing corresponding expression value (FPKM)
%
% OPTIONAL INPUTS:
%    printLevel:         Printlevel for output (default 0);
%
% OUTPUTS:
%
%   gene_id:             vector of gene identifiers present in the model
%                        that are associated with expression data
%
%   gene_expr:           vector of expression values associated to each
%                        'gened_id'
%
%   
% Authors: - S. Opdam & A. Richelle May 2017

if ~exist('printLevel','var')
    printLevel = 0;
end

%gene_id = model.genes;
gene_expr=[];
%gene_id={};
%gene_first = model.genes;
%for i=1:length(gene_first)
%    gene_id(i)=str2num(gene_first{i});
%    gene_first{i}
%    gene_id{i}=gene_first{i}
%end
gene_id=model.genes;

%size(exprData.value)
%size(exprData.gene)
%size(gene_id)
%pause
for i = 1:numel(gene_id)
%    i
    cur_ID = gene_id{i};
%    gene_id{i}
	%dataID=find(ismember(exprData.gene,cur_ID)==1);
%    find(strcmp(exprData.gene,cur_ID)==1)
    dataID=find(strcmp(exprData.gene,cur_ID)==1);
    
%    dataID
    %dataID=find(exprData.gene==cur_ID);
    %dataID
	if isempty(dataID)
    	gene_expr(i)=-1;
    elseif length(dataID)==1
    	gene_expr(i)=exprData.value(dataID);
    elseif length(dataID)>1    	
        if printLevel > 0
            disp(['Double for ',num2str(cur_ID)])
        end
    	gene_expr(i)=mean(exprData.value(dataID));
    end
end
           
end
