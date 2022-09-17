function [depGenes,EQmatrix,IndGenes] = getGeneDepMatrix(GM_matrix)
%getGeneDepMatrix
%
% Function that assesses the linear dependencies between genes in a
% metabolites-genes binary matrix.
%
%   GM_matrix (double) Boolean matrix representing the relationships 
%              metabolites and genes
%
%   EQmatrix  (logical) Matrix indicating equality relationships between
%             genes (taken from GM_matrix). Equal genes are those that are 
%             related to exactly the same metabolites in the network.
%   IndGenes  (logical) Vector indicating those genes which are not equal
%             to any other, according to GM_matrix.
%
%   Last modified.  Iván Domenzain 2019-10-28 
%

%Get relevant dimensions 
[~,G] = size(GM_matrix);
depGenes = zeros(G,1);
EQmatrix = zeros(G,G);
if issparse(GM_matrix)
    GM_matrix = full(GM_matrix);
end
%Get original rank
original_rank = rank(GM_matrix);
for i=1:G
    tempMatrix      = GM_matrix;
    geneVector      = tempMatrix(:,i);
    tempMatrix(:,i) = [];
    newRank         = rank(tempMatrix);
    if newRank<original_rank
        depGenes(i) = 1;
    end
    %Get all vectors that are equal to i-th gene
    for j=1:G
        vectorJ = GM_matrix(:,j);
        if isequal(geneVector,vectorJ)
           EQmatrix(j,i) = 1;
        end
    end
end
%Get genes that are totally independent
IndGenes = sum(EQmatrix,2);
IndGenes(IndGenes==1) = 1;
IndGenes(IndGenes>1) = 0;
end
