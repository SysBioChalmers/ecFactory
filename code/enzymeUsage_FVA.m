function FVAtable = enzymeUsage_FVA(model,enzymes)
% enzymeUsage_FVA
%   Perf
%
%   Usage: FVAtable = enzymeUsage_FVA(model,enzymes)
%
%   Ivan Domenzain.     Last edited 2020-05-27

if nargin<2
    enzymes = model.enzymes;
end
%Get parsimonious protein usages
tempModel  = model;
pool_indx  = find(strcmpi(model.rxns,'prot_pool_exchange'));
prot_indxs = find(contains(model.rxnNames,'draw_prot_'));
%prot_indxs = prot_indxs(1:(end-1));
tempModel = setParam(tempModel, 'obj',pool_indx,-1);
sol       = solveLP(tempModel,1);
%initialize variables
ranges    = zeros(length(enzymes),1);
minUsgs   = zeros(length(enzymes),1);
maxUsgs   = zeros(length(enzymes),1);
enz_pUsgs = zeros(length(enzymes),1);
enzIDs    = cell(length(enzymes),1);
if ~isempty(sol.x)
    pUsgs = sol.x(prot_indxs);
    %Loop through all the provided enzymes
    for i=1:length(enzymes)
        if ~isempty(enzymes{i})
            minFlux = 0;
            maxFlux = 0;
            pUsage  = 0;
            rxnIndx = find(contains(model.rxnNames,enzymes{i}),1);
            enzIndx = find(strcmpi(model.enzymes,enzymes{i}),1);
            enzIDs(i)  = enzymes(i);
            if ~isempty(enzIndx)
                model  = setParam(model, 'obj', rxnIndx, -1);
                sol    = solveLP(model);
                pUsage = pUsgs(enzIndx);
                if ~isempty(sol.f)
                    minFlux = sol.x(rxnIndx);
                    model   = setParam(model, 'obj', rxnIndx, +1);
                    sol     = solveLP(model);
                    if ~isempty(sol.f)
                        %disp(['Ready with enzyme #' num2str(i) ' ' model.enzymes{enzIndx}])
                        maxFlux = sol.x(rxnIndx);
                    end
                end
                ranges(i)    = (maxFlux-minFlux);
                minUsgs(i)   = minFlux;
                maxUsgs(i)   = maxFlux;
                enz_pUsgs(i) = pUsage;
            end
        else
                ranges(i)    = NaN;
                minUsgs(i)   = NaN;
                maxUsgs(i)   = NaN;
                enz_pUsgs(i) = NaN;
                %enzIDs{i}    = 'empty';
        end
    end
else
    disp('Model is not feasible')
end
varNamesT = {'enzymes' 'ranges' 'minU' 'maxU' 'pU'};
FVAtable  = table(enzIDs,ranges,minUsgs,maxUsgs,enz_pUsgs,'VariableNames', varNamesT);
end