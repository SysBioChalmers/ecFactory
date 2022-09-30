%function [optStrain,remaining,step] = ecFactory_test(model,modelParam,expYield,results_folder,graphPlot)
graphPlot = false;
model = const_ecModel;
if ~isfolder('GECKO')
    git clone --quiet --depth=1 https://github.com/SysBioChalmers/GECKO.git
end
mkdir(results_folder)
current      = pwd;
%method parameters
tol  = 1E-13; %numeric tolerance for determining non-zero enzyme usages
OEF  = 2;     %overexpression factor for enzyme targets
KDF  = 0.5;   %down-regulation factor for enzyme targets
step = 0;
%Parameters for FSEOF method
Nsteps     = 16; %number of FBA steps in ecFSEOF
alphaLims  = [0.5*expYield 2*expYield]; %biomass yield limits for ecFSEOF
lowerK     = 0.5/2;
thresholds = [0.5 1.05]; %K-score thresholds for valid gene targets
delLimit   = 0.05; %K-score limit for considering a target as deletion
%read file with essential genes list
essential = readtable('../data/essential_genes.txt','Delimiter','\t');
essential = strtrim(essential.Ids);
%Get relevant rxn indexes
modelParam.targetIndx  = find(strcmpi(model.rxns,modelParam.rxnTarget));
modelParam.CUR_indx    = find(strcmpi(model.rxns,modelParam.CSrxn));
modelParam.prot_indx   = find(strcmpi(model.rxns,'prot_pool_exchange'));
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

%2.- Add flux leak targets (those genes not optimal for production that may
%consume the product of interest. (probaly extend the approach to inmediate
%precurssors)
step = step+1;
cd (current)
disp([num2str(step) '.-  **** Find flux leak targets to block ****'])
candidates = find_flux_leaks(candidates,modelParam.targetIndx,model);
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
% 3.- discard essential genes from deletion targets
step = step+1;
disp([num2str(step) '.-  **** Removing essential genes from KD and KO targets list ****'])
[~,iB]    = ismember(candidates.genes,essential);
toRemove  = iB & candidates.k_scores<=delLimit;
candidates(toRemove,:) = [];
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')
writetable(candidates,[results_folder '/candidates_L1.txt'],'Delimiter','\t','QuoteStrings',false);
proteins = strcat('draw_prot_',candidates.enzymes);
[~,enz_pos] = ismember(proteins,model.rxns);
candidates.enz_pos = enz_pos;
% 4.- Construct Genes-metabolites network for classification of targets
step = step+1;
disp([num2str(step) '.-  **** Construct Genes-metabolites network for classification of targets ****'])
disp('  ')
%Get Genes-metabolites network
disp('  Constructing genes-metabolites graph')
disp('  ')
[GeneMetMatrix,~,Gconect] = getMetGeneMatrix(model,candidates.genes);
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

% 5.- Enzyme usage variability analysis (EVA) and prioritization of targets
step = step+1;
disp([num2str(step) '.-  **** Running EUVA for optimal production conditions (minimal biomass) ****'])
tempModel = model;
disp(' ')
disp('  - Fixed unit glucose uptake rate')
%Fix unit C source uptake
tempModel.lb(modelParam.CUR_indx) = (1-tol)*1;
tempModel.ub(modelParam.CUR_indx) = (1+tol)*1;
disp('  - Fixed suboptimal biomass production, according to provided experimental yield')
%Fix suboptimal experimental biomass yield conditions
V_bio = expYield*modelParam.CS_MW;
tempModel.lb(modelParam.growth_indx) = V_bio;
disp(['    V_bio = ' num2str(V_bio) ' h-1'])
disp('  - Production rate constrained to its maximum attainable value')
%Get and fix optimal production rate
tempModel = setParam(tempModel, 'obj', modelParam.targetIndx, +1);
sol       = solveLP(tempModel,1);
max_prod   = sol.x(modelParam.targetIndx);
tempModel.lb(modelParam.targetIndx) = (1-tol)*max_prod;
tempModel.ub(modelParam.targetIndx) = (1+tol)*max_prod;
disp(['    V_prod_max = ' num2str(max_prod) ' mmol/gDwh'])
modelParam.WT_prod = max_prod;
modelParam.V_bio   = V_bio;

disp(' ')
%Run FVA for all enzyme usages subject to fixed CUR and Grates
EVAtable = enzymeUsage_FVA(tempModel,candidates.enzymes);

disp('* Discard targets according to EUVA results for optimal production ****')
%Classify targets according to enzyme variability ranges
candidateUsages = EVAtable.pU;
disp('  ')
candidates.EV_type = cell(height(candidates),1);
candidates.EV_type(:) = {''};
idxs = find(EVAtable.minU<=tol & EVAtable.maxU<=tol);
candidates.EV_type(idxs) = {'unusable'};
idxs = find(EVAtable.minU>tol);
candidates.EV_type(idxs) = {'essential'}; %enzymes needed for optimal production
idxs = find(EVAtable.minU<=tol & EVAtable.maxU>tol);
candidates.EV_type(idxs) = {'totally_variable'}; %enzymes that can take any flux (ussually isoenzymes)
idxs = find((EVAtable.minU./EVAtable.maxU)>=0.99 & EVAtable.minU>tol);
candidates.EV_type(idxs) = {'essential_tightly_const'}; %enzymes that are needed up to its maximum capacity
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages>tol);
candidates.EV_type(idxs) = {'production_opt'};  %isoenzymes that are chosen for optimal production 
idxs = find(strcmpi(candidates.EV_type,'production_opt') & (candidateUsages./EVAtable.maxU)>=0.99);
candidates.EV_type(idxs) = {'production_opt_tight'}; %isoenzymes that are chosen for optimal production, up to its maximum capacity
idxs = find(strcmpi(candidates.EV_type,'totally_variable') & candidateUsages<=tol);
candidates.EV_type(idxs) = {'unnecessary_prod'}; %isoenzymes that are not chosen for optimal production
%Append EUVA results to target candidates table
candidates.OE(strcmpi(candidates.actions,'OE')) = OEF;
candidates.OE(strcmpi(candidates.actions,'KD')) = KDF;
candidates.OE(strcmpi(candidates.actions,'KO')) = 0;
candidates.minUsage = EVAtable.minU;
candidates.maxUsage = EVAtable.maxU;
candidates.pUsage   = candidateUsages;
%Discard enzymes 
disp('  - Discard OE targets with lb=ub=0')
toRemove = (strcmpi(candidates.EV_type,'unusable') & strcmpi(candidates.actions,'OE')) | isnan(candidates.minUsage);
candidates(toRemove,:) = [];
disp('  - Discard enzymes essential for production from deletion targets')
toRemove = (strcmpi(candidates.EV_type,'essential_tightly_const') | strcmpi(candidates.EV_type,'essential')) & ...
       (candidates.k_scores<=delLimit);
candidates(toRemove,:) = [];       
disp(' - Discard isoenzyme groups for KD/KO that contain an optimal isoform (redundant groups that increase mutant complexity)')
toRemove = [];
for k=1:max(candidates.groups)
    idx = find(candidates.groups==k & ~strcmpi(candidates.actions,'OE'));
    if length(idx)>1
        if ~isempty(candidates.enzymes(idx(1))) & any(candidates.pUsage(idx)>0)
            toRemove = [toRemove;idx];
        end
    end
end
candidates(toRemove,:) = [];
disp('  ')
disp([' * ' num2str(height(candidates)) ' gene targets remain']) 
disp('  ')

%7.- EUVA for suboptimal biomasss production subject to a minimal (1%)
% production rate of the target product and a unit CS uptake rate
%Get max biomass 
step = step+1;
disp('  ')
disp([num2str(step) '.-  **** Running EUVA for reference strain  ****'])
disp('  - Fixed unit glucose uptake rate')
disp('  - Production rate subject to a LB of 1% of its max. value')
disp('  - Biomass production fixed to its maximum attainable value')
tempModel = setParam(tempModel, 'obj', modelParam.growth_indx, +1);
tempModel.lb(modelParam.targetIndx) = 0.01*max_prod;
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
disp(' ')
%discard something 
disp('* Discard targets according to EUVA results for reference strain ****')
candidates = compare_EUVR(candidates);
disp('  - Discarding enzymes with inconsistent enzyme usage variability patterns')
toRemove = strcmpi(candidates.EUV_comparison,'embedded') | ...
           (~strcmpi(candidates.actions,'OE') & contains(candidates.EUV_comparison,'up_')) | ...
           (strcmpi(candidates.actions,'OE') & contains(candidates.EUV_comparison,'down_'));% |...
%disp(candidates(toRemove,:))
candidates(toRemove,:) = [];
disp(' ')
disp([' ' num2str(height(candidates)) ' gene targets remain']) 
disp(' ')

step = step+1;
disp([num2str(step) '.-  **** Rank targets by priority levels ****'])
disp(' ')
%Rank targets by priority
candidates.priority = zeros(height(candidates),1)+3;
candidates.priority(contains(candidates.EUV_comparison,'up_') | contains(candidates.EUV_comparison,'down_')) = 2; % second priority for genes with overlaped demand ranges
candidates.priority(contains(candidates.EUV_comparison,'_distinct')) = 1; %higher priority to the clearly up-down regulated genes
disp(['  - ' num2str(sum(candidates.priority==1)) ' targets with priority level 1 (distinct enzyme demand levels between optimal production and optimal biomass cases)']) 
disp(['  - ' num2str(sum(candidates.priority==2)) ' targets with priority level 2 (overlapped enzyme demand levels between optimal production and optimal biomass cases)']) 
disp(['  - ' num2str(sum(candidates.priority==3)) ' targets with priority level 3 (other cases)']) 
disp(' ')
%Rename enzyme variability type for KDs and KOs according to their
%variability ranges for maximum biomass production
disp(' * Classifying targets according to their enzyme usage variability range type')
idxs = ~strcmpi(candidates.actions,'OE') & candidates.pUsage < candidates.pUsageBio & candidates.pUsageBio>0 & ~contains(candidates.EV_type,'unnecessary');
candidates.EV_type(idxs) = {'biomass_opt'};
bioRatio = V_bio/maxVBio;
ratios   = candidates.pUsage./candidates.pUsageBio;
idxs     = ratios < bioRatio+1E-9 & ratios > bioRatio-1E-9;
candidates.EV_type(idxs) = {'biomass_coupled'};
disp(' ')
writetable(candidates,[results_folder '/candidates_L2.txt'],'Delimiter','\t','QuoteStrings',false);
% 8.- Combine targets
step = step+1;
disp([num2str(step) '.-  **** Find an optimal combination of remaining targets ****'])
disp(' ')
%Unconstrain CUR, unconstrain product formation 
%and set a minimal biomass formation
tempModel = setParam(tempModel,'ub',modelParam.CUR_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.CUR_indx,0);
tempModel = setParam(tempModel,'ub',modelParam.growth_indx,1000);
tempModel = setParam(tempModel,'lb',modelParam.growth_indx,0.99*V_bio);
tempModel = setParam(tempModel,'ub',modelParam.targetIndx,1000);
tempModel = setParam(tempModel,'lb',modelParam.targetIndx,0);
%set Max product formation as objective function
tempModel = setParam(tempModel,'obj',modelParam.targetIndx,+1);
%constrain enzyme usages with optimal biomass formation profile (WT)
tempModel.lb(candidates.enz_pos(find(candidates.enz_pos))) = 0.99*candidates.pUsageBio(find(candidates.enz_pos));
tempModel.ub(candidates.enz_pos(find(candidates.enz_pos))) = 1.01*candidates.maxUsageBio(find(candidates.enz_pos));
%Run mechanistic validation of targets
% [optStrain,remaining] = constructMinimalMutant(tempModel,candidates,modelParam);
% %get a gene-mets graph with the remaining targets 
% if graphPlot
%     MetsIndxs    = (~contains(optStrain.metNames,'prot_'));
%     nodeMets     = optStrain.mets(MetsIndxs);
%     toKeep       = (remaining.priority>0 & remaining.conectivity<15);
%     [GeneMetMatrix,~,~] = getMetGeneMatrix(optStrain,remaining.genes);
%     tempGMmatrix = GeneMetMatrix(:,toKeep);
%     getMGgraph(tempGMmatrix,nodeMets,tempModel,'force',remaining(toKeep,:));
% end
% disp([' * The predicted optimal strain contains ' num2str(height(remaining)) ' gene modifications']) 
% disp(' ')
% 
% cd (current) 
% writetable(remaining,[results_folder '/candidates_L3.txt'],'Delimiter','\t','QuoteStrings',false);
% %Generate transporter targets file (lists a number of transport steps
% %with no enzymatic annotation that are relevant for enhancing target
% %product formation.
% origin = 'GECKO/geckomat/utilities/ecFSEOF/results/*';
% copyfile(origin,results_folder)
% step = step+1;
% disp([num2str(step) '.-  **** Find transporter reactions with no enzyme association predicted as targets by ecFSEOF ****'])
% disp(' ')
% rxnsTable     = readtable([results_folder '/rxnsResults_ecFSEOF.txt'],'Delimiter','\t');
% transpTargets = getTransportTargets(rxnsTable,tempModel);
% disp([' * ' num2str(height(transpTargets)) ' transport reaction targets were found']) 
% disp(' ')
% writetable(transpTargets,[results_folder '/transporter_targets.txt'],'Delimiter','\t','QuoteStrings',false);
% delete([results_folder '/rxnsResults_ecFSEOF.txt'])
% delete([results_folder '/genesResults_ecFSEOF.txt'])
% %end