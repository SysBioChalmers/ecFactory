function [geneGroups,groupsCol] = getGenesGroups(EQ_Gmatrix,genes)
G = length(genes);
groupsCol   = zeros(G,1);
groupNumber = 0;
geneGroups  = cell(0);
for i=1:G
    if groupsCol(i) == 0 
    	EQgenes = find(EQ_Gmatrix(:,i));
        if numel(EQgenes)>1
            groupNumber = groupNumber+1;
            geneGroups(groupNumber) = {genes(EQgenes)};
            groupsCol(EQgenes) = groupNumber;
        end
    end
end
end