function [mutantModel,success] = getMutantModel(model,modifications,base_usage,verbose)
% getMutant
%   Get an in-Silico mutant strain based on the provided modifications
%   list. Multiple and combinatorial gene deletions, overexpressions and
%   heterologous enzyme expressions are allowed. Overexpressions are
%   applied directly on enzyme levels (upper bounds) for enzyme-constrained
%   GEMs and on stoichiometric coefficients for non-enzyme related
%   reactions.
%
%   model           Original model (either GEM or ecGEM)
%   modifications   [nx3] Cell array with the desired modifications. 
%                   - Column (1) contains strings with each of the individual 
%                   gene IDs to modify. 
%                   - Column (2) indicates the modification action; 'KO' for 
%                   deletions, 'OE' for overexpressions,'KD' for knock-downs,
%                   and 'EA' for modification of enzyme activity. 
%                   - Column (3) contains the numerical expression factor, 
%                   0 for KO, between 0 and 1 for KD and higher than 1 for
%                   OE. For 'EA' cases the new Kcat value should be
%                   provided as the expression factor in (1/s)
%   base_usage      Basal usage level for the enzyme to modify in [mmol/gDw h]
%
%   mutantModel     Mutant model with new constraints
%
%   Usage: mutantModel = getMutantModel(model,modifications,base_usage,message)
%
%   Ivan Domenzain.     Last edited 2022-08-19
if nargin<4
    verbose = false;
    if nargin<3
        base_usage = [];
    end
end

mutantModel = model;
genes2mod   = modifications(:,1);
actions     = modifications(:,2);
expF        = modifications(:,3);
pool_index  = (strcmpi(model.rxnNames,'prot_pool_exchange'));

for i=1:length(genes2mod)
    gene      = genes2mod{i};
    action    = actions{i};
    expFactor = expF{i};
    gene2modIndex = find(strcmpi(mutantModel.enzGenes,gene));
    %Knock-out mutants
    if strcmpi(action,'KO')
        %RAVEN function
        mutantModel = removeGenes(mutantModel,gene,false,false,false);
    end
    %If the gene to overexpress/knock-down has an associated enzyme then
    %act upon the enzyme usage pseudoRxn bounds
    if ~isempty(gene2modIndex)
        enzyme  = mutantModel.enzymes(gene2modIndex);
%         if ~isempty(base_usage) & expFactor<1
%             massDiff   = base_usage*mutantModel.MWs(gene2modIndex);
%             mutPmass   = (1-expFactor)*massDiff;
%             Ppool      = mutantModel.ub(pool_index)+mutPmass;
%             mutantModel.ub(pool_index) = Ppool;
%         end
        %Knock-downs and overexpressions
        if strcmpi(action,'KD') | strcmpi(action,'OE')
            %find enzyme exchange or draw reaction
            enzRxn = find(contains(mutantModel.rxnNames,enzyme));
            %If the enzRxn is bounded modify the UB
            if contains(mutantModel.rxnNames(enzRxn),'exchange') || mutantModel.ub(enzRxn)<100
                %mutantModel.ub(enzRxn) = mutantModel.ub(enzRxn)*expFactor;
                mutantModel.ub(enzRxn) = base_usage*expFactor;
                %mutantModel.lb(enzRxn) = 0;
            else
                %If the enzRxn is not bounded modify enzyme usage forcing its
                %usage lb to a given value
                if ~isempty(base_usage)
                    %OEs
                    if expFactor>1
                        mutantModel.lb(enzRxn) = base_usage*expFactor;
                        mutantModel.ub(enzRxn) = 1000;
                        %KDs
                    else
                        mutantModel.ub(enzRxn) = base_usage*expFactor;
                        mutantModel.lb(enzRxn) = 0;
                    end
                end
            end
            
            %make sure that e_i lb < e_i ub
            if mutantModel.ub(enzRxn)<= mutantModel.lb(enzRxn)
               mutantModel.lb(enzRxn) = 0.99*mutantModel.ub(enzRxn);
            end
        elseif strcmpi(action,'HE')
            enzName    = ['prot_' enzyme{1}];
            enzMetIndx = find(strcmpi(mutantModel.metNames,enzName));
            enzKcats   = find(mutantModel.S(enzMetIndx,:));
            %Avoid enzyme usage reaction
            enzKcats   = enzKcats(1:end-1);
            %Substitute kinetic  coefficients
            mutantModel.S(enzMetIndx,enzKcats) = -1/(3600*expFactor);%mutantModel.S(enzMetIndx,enzKcats)./expFactor;
        end
    %If not, modify all the stoichiometric coefficient for all
    %metabolites in the rxns encoded by the gene (unorthodox)
    else
        if strcmpi(action,'KD') | strcmpi(action,'OE')
            geneRxns = find(contains(mutantModel.grRules,gene));
            if ~isempty(geneRxns)
                mutantModel.S(:,geneRxns) = mutantModel.S(:,geneRxns).*expFactor;
            end
        end
    end
end
sol = solveLP(mutantModel);
if ~isempty(sol.x)
    success = true;
    %disp(sol.f)
    message = 'Mutant strain successfully constructed';
else
    success = false;
    message = ['WARNING: ' gene{1}  ' ' action{1} ' is not viable'];
end

if verbose 
    disp(message)
end
end