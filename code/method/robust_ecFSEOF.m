function [mutantStrain,filtered,step] = robust_ecFSEOF(model,rxnTarget,c_source,expYield,CS_MW,resultsFolder)
mkdir(resultsFolder)
current      = pwd;
tol          = 1E-12;
OE           = 2;
thresholds   = [0.5 1];
delLimit     = 0.05;
step         = 0;
essential = readtable('../../data/Essential_ORFs.txt','Delimiter','\t');
essential = strtrim(essential.Ids);

cd GECKO
%Get model parameters
cd geckomat
parameters = getModelParameters;
bioRXN     = parameters.bioRxn;
%Parameters for FSEOF method
Nsteps     = 16;
alphaLims  = [0.5*expYield 2*expYield];
%output files for genes and rxn targets
file1   = 'results/genesResults_ecFSEOF.txt';
file2   = 'results/rxnsResults_ecFSEOF.txt';
% 1.- Run FSEOF to find gene candidates
cd utilities/ecFSEOF
mkdir('results')
step = step+1;
disp([num2str(step) '.-  **** Running ecFSEOF method (from GECKO utilities) ****'])
results = run_ecFSEOF(model,rxnTarget,c_source,alphaLims,Nsteps,file1,file2);
genes   = results.geneTable(:,1);
disp('  ')
disp(['ecFSEOF yielded ' num2str(length(genes)) ' targets'])
disp('  ')
%Format results table
geneShorts = results.geneTable(:,2);
k_scores   = cell2mat(results.geneTable(:,3));
actions    = k_scores;
actions(actions<thresholds(1)) = 0;
actions(actions>1) = 1;
MWeigths           = [];
%Identify candidate genes in model enzymes
disp(' Extracting enzymatic information for target genes')
[~,iB]     = ismember(genes,model.enzGenes);
candidates = {};
pathways   = {};
for i=1:numel(iB)
    if iB(i)>0
        candidates = [candidates; model.enzymes(iB(i))];
        MWeigths   = [MWeigths; model.MWs(iB(i))];
        pathways   = [pathways; model.pathways(iB(i))];
    else
        candidates = [candidates; {''}];
        MWeigths   = [MWeigths; nan];
        pathways   = [pathways; {''}]; 
    end
end
disp('  ')
%Get results files structures
candidates = table(genes,candidates,geneShorts,MWeigths,pathways,actions,k_scores,'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'pathways' 'actions' 'k_scores'});
% Keep top results
disp(['Removing targets ' num2str(thresholds(1)) ' < K_score < ' num2str(thresholds(2))])
toKeep     = find((candidates.k_scores>=thresholds(2)|candidates.k_scores<=thresholds(1)));
candidates = candidates(toKeep,:);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')

% 2.- discard essential genes from deletion targets
step = step+1;
disp([num2str(step) '.-  **** Removing essential genes from KD and KO targets list ****'])
[~,iB]    = ismember(candidates.genes,essential);
toRemove  = iB & candidates.k_scores<=0.05;
candidates(toRemove,:) = [];
cd (current)
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[resultsFolder '/candidates_ecFSEOF.txt'],'Delimiter','\t','QuoteStrings',false);

% 3.- Enzyme usage variability analysis (EVA) and prioritization of targets

step = step+1;
disp([num2str(step) '.-  **** Running enzyme usage variability analysis ****'])
tempModel = model;
%Get relevant rxn indexes
targetIndx  = find(strcmpi(tempModel.rxns,rxnTarget));
CUR_indx    = find(strcmpi(tempModel.rxns,c_source));
growth_indx = find(strcmpi(tempModel.rxns,bioRXN));
prot_indx = find(contains(tempModel.rxns,'prot_pool'));
%Fix suboptimal experimental biomass yield conditions
V_bio = expYield*CS_MW;
tempModel.lb(growth_indx) = V_bio;
%Fix unit C source uptake
tempModel.lb(CUR_indx) = (1-tol)*1;
tempModel.ub(CUR_indx) = (1+tol)*1;
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', targetIndx, +1);
sol       = solveLP(tempModel,1);
WT_prod   = sol.x(targetIndx);
WT_CUR    = sol.x(CUR_indx);
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = (1+tol)*WT_prod;
%Run FVA for all enzyme usages subject to fixed CUR and Grates
FVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);
%sort results
candidateUsages = FVAtable.pU;
%Classify enzyme variability ranges
disp('  ')
disp(' Classifying enzyme usage variability ranges')
candidates.EV_type = cell(height(candidates),1);
idxs = find(FVAtable.minU<=tol & FVAtable.maxU<=tol);
candidates.EV_type(idxs) = {'unusable'};
idxs = find(FVAtable.minU>tol);
candidates.EV_type(idxs) = {'essential'};
idxs = find(FVAtable.minU<=tol & FVAtable.maxU>tol);
candidates.EV_type(idxs) = {'totally_variable'};
idxs = find((FVAtable.minU./FVAtable.maxU)>=0.99 & FVAtable.minU>tol);
candidates.EV_type(idxs) = {'tightly_const'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages>tol);
candidates.EV_type(idxs) = {'opt_isoenz'};
idxs = find(strcmpi(candidates.EV_type,'opt_isoenz') & (candidateUsages./FVAtable.maxU)>=0.99);
candidates.EV_type(idxs) = {'opt_isoenz_tight'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages<=tol);
candidates.EV_type(idxs) = {'subOpt_isoenz'};
%Classify overexpression types
for i=1:length(candidates.enzymes)
    if FVAtable.maxU(i)~=0 && candidates.actions(i)>0
        candidates.actions(i) = 1;
        %Enzymes that are more tightly constrained are classified as candidates
        %for overexpression by modification on Kcats
        if FVAtable.maxU(i)< OE*candidateUsages(i)
            candidates.actions(i) = 2;
        end 
    end
end
candidates.OE(candidates.actions>0)  = OE;
candidates.OE(candidates.actions==0) = 0;
candidates.minUsage = FVAtable.minU;
candidates.maxUsage = FVAtable.maxU;
candidates.pUsage   = candidateUsages;
%Discard enzymes whose usage LB = UB = 0
disp('  ')
disp(' Discard OE targets with lb=ub=0')
toRemove = strcmpi(candidates.EV_type,'unusable') & candidates.actions>0;
candidates(toRemove,:) = [];
disp(' Discard essential enzymes from deletion targets')
toRemove = (strcmpi(candidates.EV_type,'tightly_constrained') | strcmpi(candidates.EV_type,'essential')) & ...
       (candidates.k_scores<=delLimit);
candidates(toRemove,:) = [];       
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
%Generate table with FVA results
writetable(candidates,[resultsFolder '/candidates_enzUsageFVA.txt'],'Delimiter','\t','QuoteStrings',false);

% Assess genes redundancy
step = step+1;
disp([num2str(step) '.-  **** Ranking gene targets by priority level:'])
%disp([num2str(step) '.-  **** Assess genes redundancy ****'])
disp('  ')
%Get Genes-metabolites network
disp('  Constructing genes-metabolites graph')
disp('  ')
[GeneMetMatrix,~,Gconect] = getMetGeneMatrix(tempModel,candidates.genes);
%Get independent genes from GeneMetMatrix
disp('  Obtain redundant vectors in genes-metabolites graph (redundant targets)')
disp('  ')
[indGenes,G2Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%Append algebraic results to candidates table
candidates.unique = indGenes;
candidates.conectivity = Gconect.mets_number;
[~,groups]        = getGenesGroups(G2Gmatrix,candidates.genes);
candidates.groups = groups;
% Rank candidates by priority
disp(' Priority levels:')
candidates.priority = zeros(height(candidates),1);
disp('   1.- Gene candidates for OE with min. usage>0 (essential for production)')
idxs =(strcmpi(candidates.EV_type,'essential') | ...
       strcmpi(candidates.EV_type,'tightly_const')) & ...
       candidates.actions>0;
candidates.priority(idxs) = 1;
disp('   1.- Gene candidates for KO/KD with maxUsage=0 (unusable)')
idxs =strcmpi(candidates.EV_type,'unusable'); %& candidates.unique;
candidates.priority(idxs) = 1;
disp('   2.- Isoenzyme candidates for OE with pUsage=maxUsage>0 (optimal isoforms const.)')
idxs =strcmpi(candidates.EV_type,'opt_isoenz_tight') & candidates.actions>0;
candidates.priority(idxs) = 2;
disp('   2.- Suboptimal isoenzymes candidates for KO/KD (maxUsage>0 pUsage=0)')
idxs =strcmpi(candidates.EV_type,'subOpt_isoenz') & candidates.actions==0;
candidates.priority(idxs) = 2;
disp('   3.- Isoenzyme candidates for OE with pUsage>0 & pUsage<maxUsage (optimal isoforms)')
idxs =strcmpi(candidates.EV_type,'opt_isoenz') & candidates.actions>0;
candidates.priority(idxs) = 3;
disp('   3.- Groups of remaining isoenzymes that are not used for optimal production')
for k=1:max(candidates.groups)
    idx = find(candidates.groups==k & candidates.actions==0);
    if length(idx)>1
        if ~isempty(candidates.enzymes(idx(1))) & ~any(candidates.pUsage(idx)>tol)
            candidates.priority(idx) = 3;
        end
    end
end
disp('  ')
%Keep priority genes and sort them accordingly
disp(' Discard genes with priority level = 0')
candidates = candidates(candidates.priority>0,:);
candidates = sortrows(candidates,'priority','ascend');
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
writetable(candidates,[resultsFolder '/candidates_priority.txt'],'Delimiter','\t','QuoteStrings',false);
% 5.- Add flux leak targets
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Find flux leak targets to block ****'])
candidates = find_flux_leaks(candidates,targetIndx,tempModel);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 

% 6.- Mechanistic validations of FSEOF results
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Mechanistic validation of results ****'])
%Relevant rxn indexes
relIndexes = [CUR_indx, targetIndx, growth_indx];
%relax target rxn bounds
tempModel.lb(targetIndx) = (1-tol)*WT_prod;
tempModel.ub(targetIndx) = 1000;
%Unconstrain CUR and biomass formation
tempModel = setParam(tempModel,'ub',CUR_indx,1000);
tempModel = setParam(tempModel,'lb',CUR_indx,0);
tempModel = setParam(tempModel,'ub',growth_indx,1000);
tempModel = setParam(tempModel,'lb',growth_indx,0);
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',targetIndx,+1);
%Run mechanistic validation of targets
[FCs_y,FCs_p,validated]  = testAllmutants(candidates,tempModel,relIndexes);
%Discard genes with a negative impact on production yield
candidates.foldChange_yield = FCs_y; 
candidates.foldChange_pRate = FCs_p; 
candidates = candidates(validated,:);
disp(' Discard gene modifications with a negative impact on product yield or rate')
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[resultsFolder '/candidates_mech_validated.txt'],'Delimiter','\t','QuoteStrings',false);

% 7.- Find compatible combinations
step = step+1;
disp([num2str(step) '.-  **** Find compatible combinations ****'])
disp('  ')
% get optimal strain according to priority candidates
disp(' Constructing optimal strain')
disp(' ')
[mutantStrain,filtered] = getOptimalStrain(tempModel,candidates,[CUR_indx targetIndx growth_indx prot_indx],false);
cd (current) 
if ~isempty(filtered) & istable(filtered)
    disp(' ')
    step = step+1;
    disp([num2str(step) '.-  **** Discard redundant deletion targets ****'])
    disp(' ')
    %filtered = discardRedundancies(tempModel,filtered);
    actions  = cell(height(filtered),1);
    actions(filtered.actions==0 & filtered.k_scores<=delLimit)= {'KO'};
    actions(filtered.actions==0 & filtered.k_scores>delLimit)= {'KD'};
    actions(filtered.actions>0) = {'OE'};
    filtered.actions = actions;
    disp([' * ' num2str(height(filtered)) ' gene targets remain'])
    disp('  ')
    %Write final results
    writetable(filtered,[resultsFolder '/compatible_genes_results.txt'],'Delimiter','\t','QuoteStrings',false);
end
origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
copyfile(origin,resultsFolder)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FChanges_y,FChanges_p,validated] = testAllmutants(candidates,tempModel,indexes,tol)
if nargin<4
    tol = 1E-9;
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

