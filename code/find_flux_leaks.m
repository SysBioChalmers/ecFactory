function candidates = block_leaks(candidates,targetIdx,model)
%assuming that all models have an exchange reaction as objective
met    = model.metNames(find(model.S(:,targetIdx)));
cytIdx = find(strcmpi(model.compNames,'cytoplasm'));
exIdx  = find(strcmpi(model.compNames,'extracellular'));

disp([' Target molecule: ' met{1}])
metIdx = find(strcmpi(model.metNames,met) & model.metComps==cytIdx);
metExIdx = find(strcmpi(model.metNames,met) & model.metComps==exIdx);
rxnExIdxs = find(model.S(metExIdx,:)>0);

%Find reactions that consume the product in the cytoplasm
rxnIdxs  = find(model.S(metIdx,:)<0);
rxnIdxs  = setdiff(rxnIdxs,rxnExIdxs);
grRules  = model.grRules(rxnIdxs);
grRules = grRules(~cellfun(@isempty,grRules));
for i=1:length(grRules)
    grRule   = strsplit(grRules{i},' or ');
    for isoenzyme = grRule
        isoenzyme = strrep(isoenzyme,'(','');
        isoenzyme = strrep(isoenzyme,')','');
        isoenzyme = strrep(isoenzyme,' ','');
        genes = strsplit(isoenzyme{1},'and');
        presence = ismember(genes,candidates.genes);
        newTargets = genes(~presence);
        for j=1:length(newTargets)
            gene = strtrim(newTargets{j});
            enzGeneIdx = find(strcmp(model.enzGenes,gene));
            EV_type = {'flux_leak'};
            if ~isempty(enzGeneIdx)
                enzyme = model.enzymes(enzGeneIdx);
                MW = model.MWs(enzGeneIdx);
                shortName = model.enzNames(enzGeneIdx);
                pathway = model.pathways(enzGeneIdx);
            else
                enzyme = {''};
                MW = NaN;
                shortName = {''};
                pathway = {''};
            end
            newRow = [newTargets(j),enzyme,shortName,MW,pathway,'KO',0];
            candidates = [candidates;newRow];
        end
end

end

