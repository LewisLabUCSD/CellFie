function [expressionCol gene_used] = selectGeneFromGPR_modif(model, gene_names, gene_exp, parsedGPR)
% Map gene expression to reaction expression using the GPR rules. An AND
% will be replaced by MIN and an OR will be replaced by MAX.
%
% USAGE:
%   expressionCol = selectGeneFromGPR(model, gene_names, gene_exp, parsedGPR)
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
% OUTPUTS:
%   expressionCol:  reaction expression, corresponding to model.rxns.
%                   No gene-expression data and orphan reactions will
%                   be given a value of -1.
%
% AUTHOR: Anne Richelle, May 2017
    gene_used={};
    for i=1:length(model.rxns)
        gene_used{i}='';
    end
    expressionCol = -1*ones(length(model.rxns),1); %-1 means unknown/no data
    for i = 1:length(model.rxns)

        curExprArr=parsedGPR{i};
        curExpr= [];
        gene_potential=[];
        
        for j=1:length(curExprArr)
            if length(curExprArr{j})>=1
                geneID = find(ismember(gene_names,curExprArr{j}{1}));
                %gene_names
                if ~isempty(geneID) %if the gene is measured
                    [minGenevalue, minID]=min(gene_exp(geneID));
                     curExpr= [curExpr, minGenevalue]; %If there is data for any gene in 'AND' rule, take the minimum value
                    gene_potential=[gene_potential, gene_names(geneID(minID))];
                end
            end
        end
        if ~isempty(curExpr)
            [expressionCol(i), ID_max]=max(curExpr);%If there is data for any gene in the 'OR' rule, take the maximum value
            gene_used{i}=gene_potential{ID_max};
        end
    end
gene_used=gene_used;
end
