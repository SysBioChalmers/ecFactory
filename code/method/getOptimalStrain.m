function [optStrain,optGenes,FChanges,iB] = getOptimalStrain(model,candidates,rxnIndxs,protFlag)

tolerance = 0.001;

GURindex    = rxnIndxs(1);
targetIndex = rxnIndxs(2);
growth_indx = rxnIndxs(3);
prot_index  = rxnIndxs(4);
tempModel = setParam(model,'lb',targetIndex,0);
tempModel = setParam(tempModel,'obj',growth_indx,1);
%Get WT solutions
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',prot_index);
maxGrowth = WTsol(growth_indx);
tempModel = setParam(tempModel,'lb',growth_indx,0.1);
tempModel = setParam(tempModel,'obj',targetIndex,1);
%Get WT yield
[WTsol,flag] = solveECmodel(tempModel,tempModel,'pFBA',prot_index);
if flag
    if protFlag
        minIndex = prot_index;
    else
        minIndex = GURindex;
    end
    WTyield = WTsol(targetIndex)/(WTsol(minIndex));
    %Create mutants iteratively
    optStrain  = tempModel;
    FChanges   = [];
    genesFC    = [];
    counter    = 0;
    previousFC = 1;
    for i=[1 2 3]
        levelCandidates = candidates(candidates.priority==i,:);
        levelCandidates.foldChange = mean([levelCandidates.foldChange_yield,levelCandidates.foldChange_pRate],2);
        levelCandidates = sortrows(levelCandidates,'foldChange','descend');
        for j=1:length(levelCandidates.genes)
            counter = counter+1;
            gene   = levelCandidates.genes{j};
            enzyme = levelCandidates.enzymes{j};
            enzRxn = find(contains(optStrain.rxnNames,['prot_' enzyme]),1);
            if ~isempty(enzRxn)
                iterationModel = optStrain;
                iterationModel = setParam(iterationModel,'obj',enzRxn,1);
                tempSolution   = solveLP(iterationModel);
                fluxE  = tempSolution.x(enzRxn)>1E-18;
                short  = levelCandidates.shortNames{j};
                action = levelCandidates.actions(j);
                maxUse = levelCandidates.maxUsage(j);
                OEf    = levelCandidates.OE(j);
                %Avoid including enzymes that cannot carry any flux
                if fluxE & maxUse>=0
                    modifications = {gene action OEf};
                    enzUsage      = WTsol(enzRxn);
                    if action ~= 0
                        action = 2;%sol(enzRxn);%candidates.maxUsage(i);
                        numTol = tolerance;
                    else
                        numTol = 1E-12;
                    end
                    
                    if enzUsage==0
                        enzUsage = 1E-8;
                    end
                    tempMutant = getMutant(optStrain,modifications,enzUsage);
                    [mutSol,~] = solveECmodel(tempMutant,model,'pFBA',prot_index);
                    
                    if ~isempty(mutSol)
                        yield = mutSol(targetIndex)/(mutSol(minIndex));
                        FC_y  = yield/WTyield;
                        FC_p  = mutSol(targetIndex)/WTsol(targetIndex);
                        score = mean([FC_y,FC_p]);
                        %Just keep those genes that don't affect the production phenotype
                        if score >= (previousFC+numTol)
                            FChanges   = [FChanges; score];
                            genesFC    = [genesFC;{gene}];
                            optStrain  = tempMutant;
                            previousFC = score;
                            %counter = counter+1;
                            disp(['  Added target #' num2str(counter) ': (' short ')' '|  FC:' num2str(score)])
                        end
                    end
                end
            end
        end
    end
    if ~isempty(genesFC)
        [~,iB]   = ismember(genesFC,candidates.genes);
        iB       = sort(iB,'ascend');
        FChanges = table(genesFC,FChanges,'VariableNames',{'genes' 'FC'});
        optGenes = candidates(iB,:);
    else
        FChanges = table();
        optGenes = cell(1,1);
    end
end
end