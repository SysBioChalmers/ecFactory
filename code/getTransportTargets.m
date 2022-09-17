function candidates = getTransportTargets(rxnsTable,model)
candidates = table();
%identify transport and exchange rections by name
idxs = find(contains(rxnsTable.rxnNames,' transport') | contains(rxnsTable.rxnNames,' exchange'));
candidates = rxnsTable(idxs,:);
%remove found reactions from rxnsTable
rxnsTable(idxs,:) = [];
%Search for missing transport reaction by parsing rxns and identifying
%those that involve same metabolites as products and substrates
idxs = [];
for i=1:height(rxnsTable)
    rxn = rxnsTable.rxn_IDs(i);
    idx = find(strcmp(model.rxns,rxn));
    subs = model.metNames(find(model.S(:,idx)<0));
    prod = model.metNames(find(model.S(:,idx)>0));
    x = setdiff(subs,prod);
    y = setdiff(prod,subs);
    if isempty(x) & isempty(y)
        candidates = [candidates; rxnsTable(i,:)];
    end
end
end