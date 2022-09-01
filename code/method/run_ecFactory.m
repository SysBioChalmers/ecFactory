%function [mutantStrain,filtered,step] = run_ecFactory(model,modelParam,expYield,results_folder)
model = const_ecModel;
mkdir(results_folder)
current      = pwd;
%method parameters
tol  = 1E-12; %numeric tolerance for determining nin-zero enzyme usages
OEF  = 2;     %overexpression factor for enzyme targets
KDF  = 0.5;   %down-regulation factor for enzyme targets
step = 0;
%Parameters for FSEOF method
Nsteps     = 16; %number of FBA steps in ecFSEOF
alphaLims  = [0.5*expYield 2*expYield]; %biomass yield limits for ecFSEOF
thresholds = [0.5 1.05]; %K-score thresholds for valid gene targets
delLimit   = 0.05; %K-score limit for considering a target as deletion
%read file with essential genes list
essential = readtable('../../data/essential_genes.txt','Delimiter','\t');
essential = strtrim(essential.Ids);
%Get relevant rxn indexes
modelParam.targetIndx  = find(strcmpi(model.rxns,modelParam.rxnTarget));
modelParam.CUR_indx    = find(strcmpi(model.rxns,modelParam.CSrxn));
modelParam.prot_indx   = find(contains(model.rxns,'prot_pool'));
modelParam.growth_indx = find(strcmpi(model.rxns,modelParam.growthRxn));
%ecModel verification steps
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
disp(['ecFSEOF returned ' num2str(length(genes)) ' targets'])
disp('  ')
%Format results table
geneShorts = results.geneTable(:,2);
k_scores   = cell2mat(results.geneTable(:,3));
actions    = cell(numel(k_scores),1);
actions(k_scores>=thresholds(2)) = {'OE'};
actions(k_scores<thresholds(1))  = {'KD'};
actions(k_scores<delLimit) = {'KO'};
MWeigths = [];
%Identify candidate genes in model enzymes
disp(' Extracting enzymatic information for target genes')
[~,iB]     = ismember(genes,model.enzGenes);
candidates = {};
pathways   = {};
%optimize!
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
%Get results table
candidates = table(genes,candidates,geneShorts,MWeigths,pathways,actions,k_scores,'VariableNames',{'genes' 'enzymes' 'shortNames' 'MWs' 'pathways' 'actions' 'k_scores'});
%Keep results that comply with the specified K-score thresholds
disp(['Removing targets ' num2str(thresholds(1)) ' < K_score < ' num2str(thresholds(2))])
toKeep     = find((candidates.k_scores>=thresholds(2)|candidates.k_scores<=thresholds(1)));
candidates = candidates(toKeep,:);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')

% 2.- discard essential genes from deletion targets
step = step+1;
disp([num2str(step) '.-  **** Removing essential genes from KD and KO targets list ****'])
[~,iB]    = ismember(candidates.genes,essential);
toRemove  = iB & candidates.k_scores<=delLimit;
candidates(toRemove,:) = [];
cd (current)
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[results_folder '/candidates_L1.txt'],'Delimiter','\t','QuoteStrings',false);

% 3.- Enzyme usage variability analysis (EVA) and prioritization of targets
step = step+1;
disp([num2str(step) '.-  **** Running EUVA for optimal production strain (minimal biomass) ****'])
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
tempModel.lb(modelParam.targetIndx) = (1-tol)*WT_prod;
tempModel.ub(modelParam.targetIndx) = (1+tol)*WT_prod;
%Run FVA for all enzyme usages subject to fixed CUR and Grates
EVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);
%Classify targets according to enzyme variability ranges
candidateUsages = EVAtable.pU;
disp('  ')
disp(' Classifying enzyme usage variability ranges')
candidates.EV_type = cell(height(candidates),1);
idxs = find(EVAtable.minU<=tol & EVAtable.maxU<=tol);
candidates.EV_type(idxs) = {'unusable'};
idxs = find(EVAtable.minU>tol);
candidates.EV_type(idxs) = {'essential'};
idxs = find(EVAtable.minU<=tol & EVAtable.maxU>tol);
candidates.EV_type(idxs) = {'totally_variable'};
idxs = find((EVAtable.minU./EVAtable.maxU)>=0.99 & EVAtable.minU>tol);
candidates.EV_type(idxs) = {'tightly_const'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages>tol);
candidates.EV_type(idxs) = {'opt_isoenz'};
idxs = find(strcmpi(candidates.EV_type,'opt_isoenz') & (candidateUsages./EVAtable.maxU)>=0.99);
candidates.EV_type(idxs) = {'opt_isoenz_tight'};
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages<=tol);
candidates.EV_type(idxs) = {'subOpt_isoenz'};
%Append EUVA results to target candidates table
candidates.OE(strcmpi(candidates.actions,'OE')) = OEF;
candidates.OE(strcmpi(candidates.actions,'KD')) = KDF;
candidates.OE(strcmpi(candidates.actions,'KO')) = 0;
candidates.minUsage = EVAtable.minU;
candidates.maxUsage = EVAtable.maxU;
candidates.pUsage   = candidateUsages;
%Discard enzymes whose usage LB = UB = 0
disp('  ')
disp(' Discard OE targets with lb=ub=0')
toRemove = strcmpi(candidates.EV_type,'unusable') & strcmpi(candidates.actions,'OE');
candidates(toRemove,:) = [];
disp(' Discard enzymes essential for production from deletion targets')
toRemove = (strcmpi(candidates.EV_type,'tightly_constrained') | strcmpi(candidates.EV_type,'essential')) & ...
       (candidates.k_scores<=delLimit);
candidates(toRemove,:) = [];       
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
%4.- Rank targets by priority levels
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
%generate gene-gene matrix and identify linearly independent targets
[indGenes,G2Gmatrix,~] = getGeneDepMatrix(GeneMetMatrix);
%find unique targets (with no isoenzymes or not part of complexes)
candidates.unique = indGenes;
%number of metabolites connected to each gene
candidates.conectivity = Gconect.mets_number; 
%Get gene target groups (those connected to exactly the same metabolites)
[~,groups]        = getGenesGroups(G2Gmatrix,candidates.genes);
candidates.groups = groups;
% Rank candidates by priority
disp(' Priority levels:')
candidates.priority = zeros(height(candidates),1);
disp('   1.- Gene candidates for OE with min. usage>0 (essential for production)')
idxs =(strcmpi(candidates.EV_type,'essential') | ...
       strcmpi(candidates.EV_type,'tightly_const')) & ...
       strcmpi(candidates.actions,'OE');
candidates.priority(idxs) = 1;
disp('   1.- Gene candidates for KO/KD with maxUsage=0 (unusable)')
idxs =strcmpi(candidates.EV_type,'unusable'); %& candidates.unique;
candidates.priority(idxs) = 1;
disp('   2.- Isoenzyme candidates for OE with pUsage=maxUsage>0 (optimal isoforms const.)')
idxs =strcmpi(candidates.EV_type,'opt_isoenz_tight') & strcmpi(candidates.actions,'OE');
candidates.priority(idxs) = 2;
disp('   2.- Suboptimal isoenzymes candidates for KO/KD (maxUsage>0 pUsage=0)')
idxs =strcmpi(candidates.EV_type,'subOpt_isoenz') & ~strcmpi(candidates.actions,'OE');
candidates.priority(idxs) = 2;
disp('   3.- Isoenzyme candidates for OE with pUsage>0 & pUsage<maxUsage (optimal isoforms)')
idxs =strcmpi(candidates.EV_type,'opt_isoenz') & strcmpi(candidates.actions,'OE');
candidates.priority(idxs) = 3;
disp('   3.- Groups of remaining isoenzymes that are not used for optimal production')
for k=1:max(candidates.groups)
    idx = find(candidates.groups==k & ~strcmpi(candidates.actions,'OE'));
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
%5.- Add flux leak targets (those genes not optimal for production that may
%consume the product of interest. (probaly extend the approach to inmediate
%precurssors)
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Find flux leak targets to block ****'])
candidates = find_flux_leaks(candidates,modelParam.targetIndx,tempModel);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 

%6.- EUVA for suboptimal biomasss production subject to a minimal (1%)
% production rate of the target product and a unit CS uptake rate
%Get max biomass 
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Running EUVA for WT strain (max. biomass and min. product) ****'])
tempModel = setParam(tempModel, 'obj', modelParam.growth_indx, +1);
tempModel.lb(modelParam.targetIndx) = 0.01*WT_prod;
tempModel.ub(modelParam.targetIndx) = 1000;
sol       = solveLP(tempModel,1);
maxVBio   = sol.x(modelParam.growth_indx);
%Fix optimal biomass formation rate
tempModel.lb(modelParam.growth_indx) = (1-tol)*maxVBio;
tempModel.ub(modelParam.growth_indx) = (1+tol)*maxVBio;
%run EUVA for optimal biomass formation
EVAbio = enzymeUsage_FVA(tempModel,candidates.enzymes);
candidates.minUsageBio = EVAbio.minU;
candidates.maxUsageBio = EVAbio.maxU;
candidates.pUsageBio   = EVAbio.pU;

% 7.- Mechanistic validation of remaining targets
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Validation of individual targets according to FBA performance ****'])
%Unconstrain CUR, unconstrain product formation 
%and set a minimal biomass formation
tempModel = setParam(tempModel,'ub',modelParam.CUR_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.CUR_indx,0);
tempModel = setParam(tempModel,'ub',modelParam.growth_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.growth_indx,V_bio);
tempModel = setParam(tempModel,'ub',modelParam.targetIndx,1000);
tempModel = setParam(tempModel,'lb',modelParam.targetIndx,0);
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',modelParam.targetIndx,+1);
%constrain enzyme usages with optimal biomass formation profile (WT)
proteins = strcat('draw_prot_',candidates.enzymes);
[~,enz_pos] = ismember(proteins,tempModel.rxns);
candidates.enz_pos = enz_pos;
tempModel.lb(candidates.enz_pos(find(candidates.enz_pos))) = 0.99*candidates.minUsageBio(find(candidates.enz_pos));
tempModel.ub(candidates.enz_pos(find(candidates.enz_pos))) = 1.01*candidates.maxUsageBio(find(candidates.enz_pos));
%Run mechanistic validation of targets
[mutResults,topGene] = testGeneModifications(candidates,tempModel,modelParam);
%Discard genes with a negative impact on production yield or rate
candidates.performance = mutResults.performance; 
candidates = candidates(mutResults.validated,:);
disp(' Discard individual gene modifications with a negative impact on product yield or rate')
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[results_folder '/candidates_L2.txt'],'Delimiter','\t','QuoteStrings',false);

% 8.- Construct an optimal strain with remaining targets
step = step+1;
disp([num2str(step) '.-  **** Construct an optimal strain with remaining targets ****'])
disp('  ')
% get optimal strain according to priority candidates
disp(' Constructing optimal strain')
disp(' ')
[mutantStrain,filtered] = constructOptimalStrain(tempModel,candidates,modelParam);
cd (current) 
if ~isempty(filtered) & istable(filtered)
    disp(' ')
    step = step+1;
    disp([num2str(step) '.-  **** Discard redundant deletion targets ****'])
    disp(' ')
    %filtered = discardRedundancies(tempModel,filtered);
    %remove unnecessary columns    
    disp([' * ' num2str(height(filtered)) ' gene targets remain'])
    disp('  ')
    %Write final results
    writetable(filtered,[results_folder '/candidates_L3.txt'],'Delimiter','\t','QuoteStrings',false);
end
%Generate transporter targets file (lists a number of transport steps
%with no enzymatic annotation that are relevant for enhancing target
%product formation.
origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
copyfile(origin,results_folder)
rxnsTable     = readtable([results_folder '/rxnsResults_ecFSEOF.txt'],'Delimiter','\t');
transpTargets = getTransportTargets(rxnsTable,tempModel);
writetable(transpTargets,[results_folder '/transporter_targets.txt'],'Delimiter','\t','QuoteStrings',false);
delete([results_folder '/rxnsResults_ecFSEOF.txt'])
delete([results_folder '/genesResults_ecFSEOF.txt'])
%end