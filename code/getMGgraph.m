function [metGeneGraph,metNames] = getMGgraph(GeneMetMatrix,mets,model,algorithm,genes)
exclude   = {'ATP' 'ADP' 'AMP' 'GTP' 'GDP' 'GMP' 'CTP' 'CMP' 'CDP'...
             'dGDP' 'dCDP' 'dATP' 'dGTP' 'dTTP' 'dTMP' 'dGMP' 'dAMP' ...
             'dADP' 'dUMP' 'dUTP' 'diphosphate' 'phosphate' 'H+' 'H2O' ...
             'NADPH' 'NADP(+)' 'NADH' 'NAD' 'pmet_' 'carbon dioxide' 'coenzyme A' 'bicarbonate'};
tempGMmat = full(GeneMetMatrix);
[~,iB]    = ismember(mets,model.mets);
metNames  = model.metNames(iB);

if nargin<4
    algorithm = [];
end
actions = genes.actions;
kscores = genes.k_scores;
genes = genes.shortNames;
%Exclude highly connected metabolites from GeneMetMatrix
for compound=exclude
    if strcmpi(compound,'pmet_')
        metIndxs = find(contains(metNames,compound));
    else
        metIndxs = find(strcmpi(metNames,compound));
    end
    metNames(metIndxs)    = [];
    mets(metIndxs)        = [];
    tempGMmat(metIndxs,:) = [];
    iB(metIndxs) = [];
end

for i=1:length(iB)
    index = iB(i);
    comp = model.metComps(index);
    str  = ['-[' model.comps{comp} ']'];
    metNames{i} = [metNames{i} str];
end

[M,G]        = size(tempGMmat);
metGeneGraph = zeros(M,M);
for i=1:G
    geneRow  = tempGMmat(:,i);
    metIndxs = find(geneRow);
    metGeneGraph(metIndxs,metIndxs) = i;
end
%Exclude unconnected nodes from graph
toKeep = sum(logical(metGeneGraph),2)>1;
metGeneGraph = metGeneGraph(toKeep,toKeep);
metGeneGraph = graph(metGeneGraph,mets(toKeep),'OmitSelfLoops');
if nargin>4
    x = metGeneGraph.Edges.Weight;
    metGeneGraph.Edges.Name = genes(x);
end
actions = actions(x);
kscores = kscores(x);
kscores = kscores+1;
kscores = log2(kscores)+1;
if ~isempty(algorithm)
    figure                 % Creates a figure
    % Creates an axes and sets its FontSize to 18
    p = plot(metGeneGraph);
    p.NodeColor  = 'black';
    p.MarkerSize = 10;
    LWidths = kscores;

    p.LineWidth  = LWidths;
    colors = zeros(length(x),3);
    colors(strcmpi(actions,'OE'),1) = 1;
    colors(strcmpi(actions,'KD'),3) = 1;
    colors(strcmpi(actions,'KO'),3) = 0.5;  
    colors(strcmpi(actions,'KO'),2) = 0.5;
    colors(strcmpi(actions,'KO'),1) = 0.5;
    
    p.EdgeColor = colors;
    labelnode(p,1:height(metGeneGraph.Nodes),metNames(toKeep))
    layout(p,algorithm)
    labeledge(p,1:numedges(metGeneGraph),genes(x))
    title('Metabolites-genes graph')
    set(gca,'FontSize',15)
end
metNames = metNames(toKeep);
end