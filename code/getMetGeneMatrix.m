function [GeneMetMatrix,metsConectivity,genesConectivity] = getMetGeneMatrix(model,genes)
%getMetGeneMatrix
%   
%   Function that obtains a binary matrix in which rows represent
%   metabolites and columns genes. Each non-zero coeffiecient represents
%   a relationship between a gene and a metabolite (i.e. the metabolite is
%   either consumed/produced by a reaction encoded by a given gene).
%
%   model            (struct) A MATLAB GEM structure
%   genes            List of gene IDs or gene indexes to take from the model 
%
%   GeneMetMatrix    (sparse) Boolean matrix representing the relationships 
%                    metabolites and genes
%   metsConectivity  (vector) Vector that indicates the amount of
%                    metabolites related to each gene
%   genesConectivity (vector) Indicates the amount of genes related to each
%                    metabolite
%
%   Usage:  [GeneMetMatrix, metsConectivity,genesConectivity] = getGeneMetMatrix(model,genes)
% 
%   Last modified.  Iván Domenzain 2019-10-29
%

%Manage exceptions
if nargin<2
    genes = model.genes;
end

if ~isempty(genes)
    if ~isnumeric(genes)
        [iA,iB] = ismember(genes,model.genes);
        if sum(iA)~=numel(genes)
            disp('Not all provided genes were found in model')
        end
        genes = iB;
    end
else
    genes = 1:numel(model.genes);
end
%Standardize rxnGeneMat
[~,rxnGeneMat] = standardizeGrRules(model,true);
%Avoid enzyme usage reactions and pseudometabolites
MetsIndxs = find(~contains(model.metNames,'prot_'));
rxnIndxs  = find(~contains(model.rxnNames,'draw_prot_') | ~contains(model.rxns,'_REV'));
%Get metGenesMatrix
Smat  = model.S(MetsIndxs,rxnIndxs);
RGmat = rxnGeneMat(rxnIndxs,genes);
GeneMetMatrix = abs(sign(Smat))*RGmat;
GeneMetMatrix = full(GeneMetMatrix);
%Calculate row and column sums
metsConectivity  = table(model.metNames(MetsIndxs),model.mets(MetsIndxs),model.metComps(MetsIndxs),sum(logical(GeneMetMatrix),2),'VariableNames',{'metNames' 'mets' 'metComps' 'genes_number'});
genesConectivity = table(model.genes(iB),model.geneShortNames(iB),sum(logical(GeneMetMatrix),1)','VariableNames',{'genes' 'short_names' 'mets_number'});
end

