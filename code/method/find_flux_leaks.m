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
    grRule = strsplit(grRules{i},' or ');
    newgroup = max(candidates.groups)+1;
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
            unique = ((numel(grRule)==1)&(numel(genes)==1));
            
            switch unique
                case 0
                    priority = 2;
                    group = newgroup;
                case 1
                    priority = 1;
                    group = 0;
            end
            newRow = [newTargets(j),enzyme,shortName,MW,pathway,0,0,EV_type,0,0,0,0,unique,0,group,priority];
            candidates = [candidates;newRow];
        end
end

end

