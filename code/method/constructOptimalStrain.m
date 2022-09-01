function [optStrain,optGenes,FChanges,iB] = constructOptimalStrain(model,candidates,modelParam)
tolerance = 0; %numeric tolerance for performance evaluation
tempModel = model;

%Get max WT production rate
WTsol_prod = solveECmodel(tempModel,tempModel,'pFBA',modelParam.prot_indx);
%Get max WT production yield
WTsol_yield = solveECmodel(tempModel,tempModel,'pFBA',modelParam.CUR_indx);
WTprodR = WTsol_prod(modelParam.targetIndx);
WTyield = WTsol_yield(modelParam.targetIndx)/(WTsol_yield(modelParam.CUR_indx));

%Create mutants iteratively
optStrain  = tempModel;
FChanges   = [];
genesFC    = [];
counter    = 0;
previousFC = 1;
%sort targets by priority level and performance metric
%[levelCandidates,idxs] = sortrows(candidates,{'priority' 'performance'},{'ascend' 'descend'});
[levelCandidates,idxs] = sortrows(candidates,{'priority' 'performance'},{'ascend' 'descend'});
remaining = levelCandidates(idxs,:);
for j=1:length(levelCandidates.genes)
    counter = counter+1;
    gene   = levelCandidates.genes{j};
    enzyme = levelCandidates.enzymes{j};
    short  = levelCandidates.shortNames{j};
    action = levelCandidates.actions(j);
    mutF   = levelCandidates.OE(j);
    if ~isempty(enzyme)
        %pUsage = max([candidates.maxUsage(i),candidates.maxUsageBio(i)]);
        if strcmpi(levelCandidates.actions{j},'OE')
            %enzUsage = max([levelCandidates.maxUsage(j),levelCandidates.maxUsageBio(j)]);
            enzUsage = levelCandidates.maxUsageBio(j);
            if enzUsage <= 1E-15
                enzUsage = levelCandidates.maxUsage(j);
            end
        else    
            enzUsage = levelCandidates.pUsageBio(j);%min([candidates.pUsage(i),candidates.pUsageBio(i)]);
        end
    else
        if ~strcmpi(candidates.actions{j},'OE')
        	enzUsage = 1E-9;
        end
    end
    modifications = {gene action mutF};
    tempMutant = getMutantModel(optStrain,modifications,enzUsage);
    tempMutant = setParam(tempMutant,'obj',modelParam.targetIndx,1);   
%    Allow flexibilization of the other candidate enzymes
    tempMutant.ub(remaining.enz_pos(find(remaining.enz_pos))) = ...
          1.01*tempMutant.ub(remaining.enz_pos(find(remaining.enz_pos)));%levelCandidates.maxUsage(find(enz_idxs));
    tempMutant.lb(remaining.enz_pos(find(remaining.enz_pos))) = ...
          0.99*tempMutant.lb(remaining.enz_pos(find(remaining.enz_pos)));%levelCandidates.minUsage(find(enz_idxs));
    
    [mutSol_r,~] = solveECmodel(tempMutant,model,'pFBA',modelParam.prot_indx);
    [mutSol_y,~] = solveECmodel(tempMutant,model,'pFBA',modelParam.CUR_indx);
    
    if ~isempty(mutSol_r) & ~isempty(mutSol_y)
        yield = mutSol_y(modelParam.targetIndx)/(mutSol_y(modelParam.CUR_indx));
        FC_y  = yield/WTyield;
        FC_p  = mutSol_r(modelParam.targetIndx)/WTprodR;
        score = mean([FC_y,FC_p]);
        %Just keep those genes that don't affect the production phenotype
        if score >= (previousFC+tolerance)
            FChanges   = [FChanges; score];
            genesFC    = [genesFC;{gene}];
            optStrain  = tempMutant;
            previousFC = score;
            %discard enzyme of the successful candidate
            idx = find(strcmpi(remaining.genes,gene));
            remaining(idx,:) = [];
            %counter = counter+1;
            disp(['  Added target #' num2str(counter) ': (' short ')' '|  FC:' num2str(score)])
        end
    end
end
%generate output
if ~isempty(genesFC)
    [~,iB]   = ismember(genesFC,candidates.genes);
    iB       = sort(iB,'ascend');
    optGenes = candidates(iB,:);
    optGenes.FC = FChanges;
    FChanges = table(genesFC,FChanges,'VariableNames',{'genes' 'FC'});
else
    FChanges = table();
    optGenes = cell(1,1);
end
end