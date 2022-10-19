function [expressionCol, gene_used] = selectGeneFromGPR_all(model, gene_names, gene_exp, parsedGPR, minSum)
% Map gene expression to reaction expression using the GPR rules. An AND
% will be replaced by MIN and an OR will be replaced by MAX.
%
% USAGE:
%   expressionCol = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR, minMax)
%
% INPUTS:
%   model:          COBRA model struct
%   gene_names:     gene identifiers corresponding to gene_exp. Names must
%                   be in the same format as model.genes (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   gene_exp:       gene FPKM/expression values, corresponding to names (column vector)
%                   (as returned by "findUsedGeneLevels.m")
%   parsedGPR:      GPR matrix as returned by "GPRparser.m"
%
% OPTIONAL INPUTS:
%   minSum:         instead of using min and max, use min for AND and Sum
%                   for OR
%
% OUTPUTS:
%   expressionCol:  reaction expression, corresponding to model.rxns.
%                   No gene-expression data and orphan reactions will
%                   be given a value of -1.
%
% AUTHOR: Anne Richelle, May 2017


if ~exist('minSum','var')
    minSum = false;
end
numSamp = length(gene_exp(1,:));
gene_used={};
for i=length(model.rxns):-1:1
    for zz = 1:numSamp
	gene_used{i,zz}='';
    end
end
    
SampStr = {};
for zz = 1:numSamp
SampStr{zz} = int2str(zz);
end

% -1 means unknown/no data
expressionCol = -1*ones(length(model.rxns),length(gene_exp(1,:))); 
for i = 1:length(model.rxns)
    curExprArr=parsedGPR{i};
    for zz = 1:numSamp
%         tmpzz = int2str(zz);
    tmpStuc.(['curExpr',SampStr{zz}])= [];
    tmpStuc.(['gene_potential',SampStr{zz}])=[];
    end
    for j=1:length(curExprArr)
        if length(curExprArr{j})>=1
            geneID = find(ismember(gene_names,curExprArr{j}));
            %geneID = find(ismember(gene_names,str2num(curExprArr{j}{1})));
            % if the gene is measured
            if ~isempty(geneID)
                for zz = 1:numSamp
                   % tmpzz = int2str(zz);
                if minSum
                    % This is an or rule, so we sum up all options.
                    tmpStuc.(['curExpr',SampStr{zz}])= [tmpStuc.(['curExpr',SampStr{zz}]), sum(gene_exp(geneID,zz))]; 
                    tmpStuc.(['gene_potential',SampStr{zz}])=[tmpStuc.(['gene_potential',SampStr{zz}]), gene_names(geneID,zz)'];
                else
                    % If there is data for any gene in 'AND' rule, take the minimum value
                    [minGenevalue, minID]=min(gene_exp(geneID,zz));
                    tmpStuc.(['curExpr',SampStr{zz}])= [tmpStuc.(['curExpr',SampStr{zz}]), minGenevalue]; %If there is data for any gene in 'AND' rule, take the minimum value
                    tmpStuc.(['gene_potential',SampStr{zz}])=[tmpStuc.(['gene_potential',SampStr{zz}]), gene_names(geneID(minID))];
                end
                end
            end
        end
    end
    for zz = 1:numSamp
%         tmpzz = int2str(zz);
    if ~isempty(tmpStuc.(['curExpr',SampStr{zz}]))
        if minSum
            % in case of min sum these are and clauses that are combined, so its the minimum.
            [expressionCol(i,zz), ID_min]=min(tmpStuc.(['curExpr',SampStr{zz}]));
            gene_used{i,zz}=tmpStuc.(['gene_potential',SampStr{zz}])(ID_min);
        else
            % if there is data for any gene in the 'OR' rule, take the maximum value
            [expressionCol(i,zz), ID_max]=max(tmpStuc.(['curExpr',SampStr{zz}]));
            gene_used{i,zz}=tmpStuc.(['gene_potential',SampStr{zz}])(ID_max);
        end
    end
    end
end

end