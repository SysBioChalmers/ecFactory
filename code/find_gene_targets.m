%% 1.- Clone GitHub repositories and load ecModel
clear
%Probably point to specific versions (the same of the publication)
cd method
git clone --quiet --depth=1 https://github.com/SysBioChalmers/GECKO.git
cd ..
%ecModels catalogue
git clone --quiet --depth=1 https://github.com/SysBioChalmers/ecModels.git

%% 2.- load ecYeastGEM 
load(['ecModels/ecYeastGEM/model/' model_name '_batch.mat'])
ecModel = ecModel_batch;
%% 3.- Set case-specific variables
%Indicate carbon source uptake reaction name
c_source   = 'D-glucose exchange (reversible)';
%Indicate carbon source mol. weight in grams/mmol
CS_MW      = 0.18015;
%Indicate the name of the growth (or biomass exchange) pseudoreaction
growth_rxn = 'growth';
%Indicate the wild-type ecModel file name (from ecModels catalogue) 
model_name = 'ecYeastGEM';
%get current directory
current = pwd;
%product name provide a product name
product_name = '2-phenylethanol';
%results folder name
results_tag    = '2_phenylethanol';
results_folder = ['../../results/' results_tag '_targets'];
mkdir(results_folder)

%% 4.- Identify production target and verify production

%Check presence of product in ecModel.metNames
prod_pos = find(strcmpi(ecModel.metNames,product_name));
%get extracellular compartment Index
ext = find(strcmpi(ecModel.compNames,'extracellular'));
if ~isempty(prod_pos)
    %get secretion reactions
    [~,exch_rxns] = getExchangeRxns(ecModel,'out');
    target_pos    = []; i = 1;
    while isempty(target_pos) & i <= numel(prod_pos)     
        target_pos = exch_rxns(ecModel.S(prod_pos(i),exch_rxns)~=0);
        i = i +1;
    end
    
    if isempty(target_pos)
    	warning(['No exchange reaction was found for: ' product_name])
        %get the extracellular met index
        ext_pos                    = prod_pos(find(ecModel.metComps(prod_pos) == ext));
        [prod_ecModel,target_name] = addExchangeRxns(ecModel,'out',ext_pos);
        target_pos                 = find(contains(prod_ecModel.rxns,target_name{1}));
        disp(['* Production reaction created for ' product_name])
    else    
    	disp(['* Production reaction found: ' ecModel.rxnNames{target_pos}])
    	formula = constructEquations(ecModel,target_pos);
        disp(formula{1})
        prod_ecModel = ecModel;
    end
    
    if isempty(target_pos)
       error(['Exchange reaction for ' product_name ' could not be created'])
    else
        target_rxn = prod_ecModel.rxns(target_pos);
        prod_ecModel = setParam(prod_ecModel,'obj',target_pos,1);
        disp([ecModel.rxnNames{target_pos} ' has been set as cellular objective'])
    end
else
    error('Your provided product is not part of the used ecModel metabolites list')
end
%Open all exchanges
[~,exch_rxns] = getExchangeRxns(ecModel,'out');
prod_ecModel.ub(exch_rxns) = 1000;
[~,exch_rxns] = getExchangeRxns(prod_ecModel,'in');
prod_ecModel.ub(exch_rxns) = 1000;
%
[~,exch_rxns] = getExchangeRxns(ecModel,'out');
%Check flux through target reacion
flux = haveFlux(prod_ecModel,1-12,target_pos);
if flux 
   disp(['* Your ecModel is suited for production of: ' product_name])
else
    error('Your ecModel is not suited for production, please revise your production pathway connectivity')
end
   
%% 5.-Set constraints

%Set media conditions for batch growth
const_ecModel = changeMedia_batch(prod_ecModel,c_source,'Min'); %model-specific script
CS_index      = find(strcmpi(const_ecModel.rxnNames,c_source));
growthPos     = find(strcmpi(const_ecModel.rxnNames,growth_rxn));
%Unconstrain growth
const_ecModel = setParam(const_ecModel,'lb',growthPos,0);
const_ecModel = setParam(const_ecModel,'ub',growthPos,1000);
%Get biomass yield for a unit glucose uptake rate
const_ecModel = setParam(const_ecModel,'obj',growthPos,1);
const_ecModel = setParam(const_ecModel,'ub',CS_index,1);
solution = solveLP(const_ecModel,1);
%Return objective to target production reaction
%const_ecModel = setParam(const_ecModel,'obj',target_pos,1);
%Check if model can carry flux for the target rxn

% Get max. biomass yield
WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
%Obtain a suboptimal yield value to run ecFactory
expYield = 0.49*WT_yield;
disp(['* The ecFactory method will scan flux distributions with a fixed biomass yield. spanning from: ' num2str(0.5*expYield) ' to: ' num2str(num2str(2*expYield)) ' g biomass/g glucose]'])

% %% 6.- Run Enzyme control analysis
% cd method
% %Run Enzyme sensitivity analysis
% try
%     %Set suboptimal growth rate as lb, target rxn as objective and
%     %unconstrained glucose uptake
%     tempModel = setParam(const_ecModel,'lb',growthPos,0.5*solution.x(growthPos));
%     tempModel = setParam(tempModel,'obj',target_pos,1);
%     tempModel = setParam(tempModel,'ub',CS_index,1000);
%     ECCs      = getECCs(tempModel,target_pos);
%     writetable(ECCs,['../../results/' results_tag '_targets/' results_tag '_ECCs.txt'],'QuoteStrings',false,'Delimiter','\t')
% catch
%     disp('The model is not suitable for ECC analysis')
% end

%% 9.- Run ecFactory method
try
    const_ecModel = check_enzyme_fields(const_ecModel);
    [optStrain,candidates,step] = ecFactory(const_ecModel,target_rxn,const_ecModel.rxns(CS_index),expYield,CS_MW,results_folder);
    candidates.FCC = zeros(height(candidates),1);
    try
        [iA,iB] = ismember(candidates.enzymes,ECCs.enzymes);
        candidates.FCC(find(iA)) = ECCs.CC(iB);
        writetable(candidates,[results_folder '/compatible_genes_results.txt'],'Delimiter','\t','QuoteStrings',false);
    catch
        disp('No ECCs file available for this case')
    end
    %Generate transporter targets file (lists a number of transport steps
    %with no enzymatic annotation that are relevant for enhancing target
    %product formation.
    rxnsTable     = readtable([results_folder '/rxnsResults_ecFSEOF.txt'],'Delimiter','\t');
    transpTargets = getTransportTargets(rxnsTable,tempModel);
    writetable(transpTargets,[results_folder '/transporter_targets.txt'],'Delimiter','\t','QuoteStrings',false);
catch
    disp('The model is not suitable for robust ecFSEOF')
end
disp(' ')
cd (current)
%%
