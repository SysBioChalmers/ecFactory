function [results,TOPgene] = testGeneModifications(candidates,tempModel,modelParam,message)
if nargin<4
    message = true;
end
tolerance  = 1E-12;
FChanges_y = [];
FChanges_p = [];
%tempModel = original;

%Get max. production rate
[WTsol1,~] = solveECmodel(tempModel,tempModel,'pFBA',modelParam.prot_indx);
maxRateWT  = WTsol1(modelParam.targetIndx);
%printFluxes(tempModel,WTsol1,true)
%Get max. production yield
[WTsol2,~] = solveECmodel(tempModel,tempModel,'pFBA',modelParam.CUR_indx);
maxYieldWT = WTsol2(modelParam.targetIndx)/WTsol2(modelParam.CUR_indx);

for i=1:height(candidates)
    modifications = {candidates.genes{i} candidates.actions{i} candidates.OE(i)};
    pUsage = [];
    if ~isempty(candidates.enzymes{i})
        %pUsage = max([candidates.maxUsage(i),candidates.maxUsageBio(i)]);
        if strcmpi(candidates.actions{i},'OE')
            pUsage = candidates.maxUsageBio(i);
            if pUsage <= 1E-15
                pUsage = candidates.maxUsage(i);
            end
        else    
            pUsage = candidates.pUsageBio(i);%min([candidates.pUsage(i),candidates.pUsageBio(i)]);
        end
    else
        if ~strcmpi(candidates.actions{i},'OE')
        	pUsage = 1E-9;
        end
    end
    mutantModel = getMutantModel(tempModel,modifications,pUsage);
    %get mutant solutions
    
    %max rate
    [mut_sol_resp,~] = solveECmodel(mutantModel,mutantModel,'pFBA',modelParam.CUR_indx);
    if ~isempty(mut_sol_resp) 
        yield = mut_sol_resp(modelParam.targetIndx)/mut_sol_resp(modelParam.CUR_indx);
        FC_y  = yield/maxYieldWT;
    else
        FC_y = 0;
    end
    %max yield
    [mut_sol_ferm,~] = solveECmodel(mutantModel,mutantModel,'pFBA',modelParam.prot_indx);
    if ~isempty(mut_sol_ferm) 
        FC_p  = mut_sol_ferm(modelParam.targetIndx)/maxRateWT;
    else
        FC_p = 0;
    end
    FChanges_y = [FChanges_y; FC_y];
    FChanges_p = [FChanges_p; FC_p];
    %disp(['Ready with genetic modification #' num2str(i) '[' short ': ' num2str(action) '] FC: ' num2str(FC_y)])
end
%Compute performance metric
performance = mean([FChanges_y,FChanges_p],2);
validated   = performance>=(1-tolerance);
%performance = (performance-1)/0.01;
results     = table(candidates.shortNames,FChanges_p,FChanges_y,performance,validated,'VariableNames',{'shortNames' 'max_rate' 'max_yield' 'performance' 'validated'});
%identify top genetic modification
[topPerf,I] = max(results.performance);
TOPenz = candidates.shortNames{I};
TOPgene = candidates.genes{I};
if message 
    disp([' Top candidate gene: ' TOPenz ' performance: ' num2str(topPerf)])
end
end