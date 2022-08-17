function [FChanges_y,FChanges_p,validated] = testMutants(candidates,tempModel,indexes,tol)
if nargin<4
    tol = 0;
end
FChanges_y = [];
FChanges_p = [];
CUR_indx    = indexes(1);
targetIndx  = indexes(2);
growth_indx = indexes(3);
%Index to minimize (bi-level optimization)
minIndex  = find(contains(tempModel.rxnNames,'prot_pool'));
tempModel = setParam(tempModel,'lb',targetIndx,0);
tempModel = setParam(tempModel,'obj',growth_indx,1);
%Get WT solutions
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',minIndex);
maxGrowth = WTsol(growth_indx);
tempModel = setParam(tempModel,'lb',growth_indx,0.1);
tempModel = setParam(tempModel,'obj',targetIndx,1);
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',minIndex);
prodWT    = WTsol(targetIndx);
WTval     = prodWT/WTsol(CUR_indx);

for i=1:height(candidates)
    gene   = candidates.genes{i};
    enzyme = candidates.enzymes{i};
    short  = candidates.shortNames{i};
    action = candidates.actions(i);
    if action ==1
        action = 2;
    end
    OEf    = candidates.OE(i);
    modifications = {gene action OEf};
    if ~isempty(enzyme)
        index  = find(contains(tempModel.rxnNames,['draw_prot_' enzyme]));
        pUsage = WTsol(index);
    else 
        pUsage = [];
    end
    if action ==0 & isempty(pUsage)
        pUsage = 1E-9;
    end
    mutantModel     = getMutant(tempModel,modifications,pUsage);
    [mutSolution,~] = solveECmodel(mutantModel,mutantModel,'pFBA',minIndex);
    if ~isempty(mutSolution)
        yield = mutSolution(targetIndx)/mutSolution(CUR_indx);
        FC_y  = yield/WTval;
        FC_p  = mutSolution(targetIndx)/prodWT;
    else
        FC_y = 0;
        FC_p = 0;
    end
    FChanges_y = [FChanges_y; FC_y];
    FChanges_p = [FChanges_p; FC_p];
    %disp(['Ready with genetic modification #' num2str(i) '[' short ': ' num2str(action) '] FC: ' num2str(FC_y)])
end
validated  = mean([FChanges_y,FChanges_p],2)>=(1-tol);
[maxVal,I] = max(FChanges_y);
if ~(maxVal>=1)
    maxVal = [];
    gene   = [];
else 
    TOPgene = candidates.genes{I};
    FC_y    = FChanges_y(I);
    disp([' Top candidate gene: ' short ' FC: ' num2str(FC_y)])
end
end

