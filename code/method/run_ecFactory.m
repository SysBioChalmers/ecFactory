function [mutantStrain,filtered,step] = run_ecFactory(model,modelParam,expYield,resultsFolder)
mkdir(resultsFolder)
current      = pwd;
tol          = 1E-12;
OE           = 2;
thresholds   = [0.5 1];
delLimit     = 0.05;
step         = 0;
essential = readtable('../../data/essential_genes.txt','Delimiter','\t');
essential = strtrim(essential.Ids);


%Get relevant rxn indexes
modelParam.targetIndx  = find(strcmpi(model.rxns,modelParam.rxnTarget));
modelParam.CUR_indx    = find(strcmpi(model.rxns,modelParam.CSrxn));
modelParam.prot_indx   = find(contains(model.rxns,'prot_pool'));
modelParam.growth_indx = find(strcmpi(model.rxns,modelParam.growthRxn));

%verification steps
model = check_enzyme_fields(model);
if ~isempty(modelParam.targetIndx)
    %Check if model can carry flux for the target rxn
    flux = haveFlux(model,1-12,modelParam.rxnTarget);
    if flux
        disp(['* Your ecModel can carry flux through the reaction: ' model.rxnNames{modelParam.targetIndx}])
    else
        disp(['* Your ecModel cannot carry flux through the reaction: ' model.rxnNames{modelParam.targetIndx} ', please check the applied constraints'])
    end
else
    error('The provided target reaction is not part of the ecModel.rxns field')
end

%Parameters for FSEOF method
Nsteps     = 16;
alphaLims  = [0.5*expYield 2*expYield];
%output files for genes and rxn targets
file1   = 'results/genesResults_ecFSEOF.txt';
file2   = 'results/rxnsResults_ecFSEOF.txt';
% 1.- Run FSEOF to find gene candidates
cd GECKO/geckomat/utilities/ecFSEOF
mkdir('results')
step = step+1;
disp([num2str(step) '.-  **** Running ecFSEOF method (from GECKO utilities) ****'])
results = run_ecFSEOF(model,modelParam.rxnTarget,modelParam.CSrxn,alphaLims,Nsteps,file1,file2);
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
writetable(candidates,[resultsFolder '/candidates_L1.txt'],'Delimiter','\t','QuoteStrings',false);

% 3.- Enzyme usage variability analysis (EVA) and prioritization of targets

step = step+1;
disp([num2str(step) '.-  **** Running enzyme usage variability analysis ****'])
tempModel = model;
%Fix suboptimal experimental biomass yield conditions
V_bio = expYield*modelParam.CS_MW;
tempModel.lb(modelParam.growth_indx) = V_bio;
%Fix unit C source uptake
tempModel.lb(modelParam.CUR_indx) = (1-tol)*1;
tempModel.ub(modelParam.CUR_indx) = (1+tol)*1;
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', modelParam.targetIndx, +1);
sol       = solveLP(tempModel,1);
WT_prod   = sol.x(modelParam.targetIndx);
WT_CUR    = sol.x(modelParam.CUR_indx);
tempModel.lb(modelParam.targetIndx) = (1-tol)*WT_prod;
tempModel.ub(modelParam.targetIndx) = (1+tol)*WT_prod;
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
%writetable(candidates,[resultsFolder '/candidates_enzUsageFVA.txt'],'Delimiter','\t','QuoteStrings',false);

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
%writetable(candidates,[resultsFolder '/candidates_priority.txt'],'Delimiter','\t','QuoteStrings',false);
% 5.- Add flux leak targets
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Find flux leak targets to block ****'])
candidates = find_flux_leaks(candidates,modelParam.targetIndx,tempModel);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 

% 6.- Mechanistic validations of FSEOF results
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Mechanistic validation of results ****'])
%relax target rxn bounds
tempModel.lb(modelParam.targetIndx) = (1-tol)*WT_prod;
tempModel.ub(modelParam.targetIndx) = 1000;
%Unconstrain CUR and biomass formation
tempModel = setParam(tempModel,'ub',modelParam.CUR_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.CUR_indx,0);
tempModel = setParam(tempModel,'ub',modelParam.growth_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.growth_indx,0);
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',modelParam.targetIndx,+1);
%Run mechanistic validation of targets
[FCs_y,FCs_p,validated]  = testMutants(candidates,tempModel,modelParam);
%Discard genes with a negative impact on production yield
candidates.foldChange_yield = FCs_y; 
candidates.foldChange_pRate = FCs_p; 
candidates = candidates(validated,:);
disp(' Discard gene modifications with a negative impact on product yield or rate')
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[resultsFolder '/candidates_L2.txt'],'Delimiter','\t','QuoteStrings',false);

% 7.- Find compatible combinations
step = step+1;
disp([num2str(step) '.-  **** Find compatible combinations ****'])
disp('  ')
% get optimal strain according to priority candidates
disp(' Constructing optimal strain')
disp(' ')
[mutantStrain,filtered] = getOptimalStrain(tempModel,candidates,modelParam,false);
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
    %remove unnecessary columns
    
    disp([' * ' num2str(height(filtered)) ' gene targets remain'])
    disp('  ')
    %Write final results
    writetable(filtered,[resultsFolder '/candidates_L3.txt'],'Delimiter','\t','QuoteStrings',false);
end
%Generate transporter targets file (lists a number of transport steps
%with no enzymatic annotation that are relevant for enhancing target
%product formation.
origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
copyfile(origin,resultsFolder)
rxnsTable     = readtable([resultsFolder '/rxnsResults_ecFSEOF.txt'],'Delimiter','\t');
transpTargets = getTransportTargets(rxnsTable,tempModel);
writetable(transpTargets,[resultsFolder '/transporter_targets.txt'],'Delimiter','\t','QuoteStrings',false);
delete([resultsFolder '/rxnsResults_ecFSEOF.txt'])
delete([resultsFolder '/genesResults_ecFSEOF.txt'])
end