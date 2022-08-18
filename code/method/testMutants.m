function [FChanges_y,FChanges_p,validated] = testMutants(candidates,tempModel,modelParam,subOptGrowth)
if nargin<4
    subOptGrowth = 0.1;
end
FChanges_y = [];
FChanges_p = [];

%unconstrain production, set a suboptimal growth rate and set production as
%an objective to maximize
tempModel = setParam(tempModel,'lb',modelParam.targetIndx,0);
tempModel = setParam(tempModel,'lb',modelParam.growth_indx,subOptGrowth);
tempModel = setParam(tempModel,'obj',modelParam.targetIndx,1);
%Get WT solutions
[WTsol,~] = solveECmodel(tempModel,tempModel,'pFBA',modelParam.prot_indx);
prodWT    = WTsol(modelParam.targetIndx);
WTval     = prodWT/WTsol(modelParam.CUR_indx);

for i=1:height(candidates)
    gene   = candidates.genes{i};
    enzyme = candidates.enzymes{i};
    OEf = candidates.OE(i);
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
    [mutSolution,~] = solveECmodel(mutantModel,mutantModel,'pFBA',modelParam.prot_indx);
    if ~isempty(mutSolution)
        yield = mutSolution(modelParam.targetIndx)/mutSolution(modelParam.CUR_indx);
        FC_y  = yield/WTval;
        FC_p  = mutSolution(modelParam.targetIndx)/prodWT;
    else
        FC_y = 0;
        FC_p = 0;
    end
    FChanges_y = [FChanges_y; FC_y];
    FChanges_p = [FChanges_p; FC_p];
    %disp(['Ready with genetic modification #' num2str(i) '[' short ': ' num2str(action) '] FC: ' num2str(FC_y)])
end
validated  = mean([FChanges_y,FChanges_p],2)>=1;
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

