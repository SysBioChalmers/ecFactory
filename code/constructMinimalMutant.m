function [optMutant,remaining] = constructMinimalMutant(model,candidates,modelParam)
tol =-1E-12;%-1E-6;
tempModel = model;
%Get max WT production rate
WTsol_prod = solveECmodel(tempModel,tempModel,'pFBA',modelParam.prot_indx);
%Get max WT production yield
WTsol_yield = solveECmodel(tempModel,tempModel,'pFBA',modelParam.CUR_indx);
WTprodR = WTsol_prod(modelParam.targetIndx);
WTyield = WTsol_yield(modelParam.targetIndx)/(WTsol_yield(modelParam.CUR_indx));
%get mutant with all modifications
optMutant=model;
toRemove = [];
for i=1:height(candidates)
    gene   = candidates.genes(i);
    action = candidates.actions(i);
    mutF   = 1;%candidates.OE(i);
    if strcmpi(action,'OE')
        enzUsage = candidates.maxUsage(i);
        if enzUsage <= 1E-15
            enzUsage = candidates.maxUsage(i);
        end
    else
        enzUsage = candidates.pUsage(i);
        if strcmpi(action,'KO')
            action = 'KD';
            enzUsage = 0;
        end
    end
    modifications = {gene action mutF};

    [tempModel,success] = getMutantModel(optMutant,modifications,enzUsage);
    if success 
       optMutant = tempModel;
    else
        toRemove = [toRemove; i];
    end
end
if ~isempty(toRemove)
    disp('The following gene modifications are not compatible with the rest of remaining candidate targets')
    disp(candidates(toRemove,[1 2 3 6]))
    candidates(toRemove,:) = [];
end
optMutant = setParam(optMutant,'obj',modelParam.targetIndx,1);
%obtain optimal production rate and yield
[mutSol_r,~] = solveECmodel(optMutant,model,'pFBA',modelParam.prot_indx);
[mutSol_y,~] = solveECmodel(optMutant,model,'pFBA',modelParam.CUR_indx);
OptprodR = mutSol_r(modelParam.targetIndx);
Optyield = mutSol_y(modelParam.targetIndx)/(mutSol_y(modelParam.CUR_indx));
bYield    = mutSol_y(modelParam.growth_indx)/(mutSol_y(modelParam.CUR_indx)*modelParam.CS_MW);
disp('Finding a minimal combination of targets displaying:')
disp([' - a production rate of: ' num2str(OptprodR) ' mmol/gDwh'])
disp([' - a production yield of: ' num2str(Optyield) ' mmol/mmol glucose'])
disp([' - a biomass yield of: ' num2str(bYield) 'g biomass/g glucose'])
disp(' ')
%sort targets by priority level and k_score
[levelCandidates,~] = sortrows(candidates,{'priority' 'k_scores'},{'ascend' 'ascend'});
counter   = 0;
remaining = table();
for i=1:height(levelCandidates)
    %reverse modifications
    gene   = levelCandidates.genes(i);
    action = levelCandidates.actions(i);
    enzIdx = levelCandidates.enz_pos(i);
    short  = levelCandidates.shortNames{i};
    mutF = 1;
    tempMutant = optMutant;
    %revert mutation
    if enzIdx>0
        saturationOpt         =  mutSol_r(enzIdx)/(optMutant.ub(enzIdx)+1E-15);
        tempMutant.ub(enzIdx) = model.ub(enzIdx);
        tempMutant.lb(enzIdx) = model.lb(enzIdx);
        if tempMutant.ub(enzIdx) <=tempMutant.lb(enzIdx)
            tempMutant.lb(enzIdx) = 0.99*tempMutant.ub(enzIdx);
        end
        
    %for reactions without enzymatic reaction
    else
        enzUsage = 1E-12;
        if strcmpi(action,'KO')
            reversal = {'OE'};
        else
            if strcmpi(action,'KD')
                reversal = {'OE'};
            else
                reversal = {'KD'};
            end
        end
        modifications = {gene reversal mutF};
        tempMutant    = getMutantModel(optMutant,modifications,enzUsage);
        saturationOpt = NaN;
    end
    %Get max WT production rate
    mutsol_prod = solveECmodel(tempMutant,tempMutant,'pFBA',modelParam.prot_indx);
    %Get max WT production yield
    mutsol_yield = solveECmodel(tempMutant,tempMutant,'pFBA',modelParam.CUR_indx);
    mutprodR     = mutsol_prod(modelParam.targetIndx);
    mutyield     = mutsol_yield(modelParam.targetIndx)/(mutsol_yield(modelParam.CUR_indx));
    if enzIdx>0 && numel(enzIdx)==1
        saturationM =  mutsol_yield(enzIdx)/(tempMutant.ub(enzIdx)+1E-15);
    else 
        saturationM = NaN;
    end
    FC_y  = mutyield/Optyield;
    FC_p  = mutprodR/OptprodR;
    score = mean([FC_y,FC_p]);
    %Discard genes that don't affect the optimal phenotype
    flag = true;
    %discard OE targets that show low saturation after reversing the
    %modification
    if (strcmpi(action,'OE') & saturationM <= (model.ub(enzIdx)/levelCandidates.maxUsage(i)))
        flag = false;
    end
    
    if isnan(score)
        score = 0;
    end
    
    if flag && ...
       (...
        (score<=1+tol) || ...
        (strcmpi(action,'OE') && saturationOpt >= 0.99) ...%%|| ((~strcmpi(action,'OE') && saturationOpt <= saturationM)) ...
       )
        
        remaining = [remaining;levelCandidates(i,:)];
        counter = counter+1;
        disp(['  Validated optimal target # ' num2str(counter) ': (' short '), ' action{1}])
    else
    	optMutant = tempMutant;
    end
end
remaining = sortrows(remaining,{'priority' 'k_scores'},{'ascend' 'descend'});
remaining= removevars(remaining,{'enz_pos' 'OE' 'minUsage' 'maxUsage' 'pUsage' 'minUsageBio' 'maxUsageBio' 'pUsageBio'});
end