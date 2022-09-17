function model = check_enzyme_fields(model)
cytIdx = find(strcmpi(model.compNames,'cytoplasm'));

if isfield(model,'enzNames')
    if length(model.enzymes)>length(model.enzNames)
        i = length(model.enzNames)+1;
        for j=i:length(model.enzymes)
            %Add enzNames
            model.enzNames(j) = model.enzGenes(j);
            %check if the draw rxn has been named properly
            enzyme = ['prot_' model.enzymes{j}];
            metIdx = find(strcmpi(model.mets,enzyme));
            rxnIdx = find(model.S(metIdx,:)>0);
            rxnName = ['draw_prot_' enzyme];
            model.rxns{rxnIdx} = rxnName;
            model.rxnNames{rxnIdx} = rxnName;   
            %Assign new enzymes to cytoplasm
            model.metComps(metIdx) = cytIdx;
            %Add missing molecular weights with average MW
            if length(model.enzymes)>length(model.MWs)
                model.MWs(j) = mean(model.MWs);
            end
             if length(model.enzymes)>length(model.pathways)
                model.pathways{j} = '';
            end
        end
    end
end
end