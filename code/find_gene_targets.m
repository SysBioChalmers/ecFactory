%% 1.- Clone GitHub repositories and load ecModel
clear
%Probably point to specific versions (the same of the publication)
cd method
git clone --quiet --depth=1 https://github.com/SysBioChalmers/GECKO.git
cd ..
%ecModels catalogue
git clone --quiet --depth=1 https://github.com/SysBioChalmers/ecModels.git
%% 2.- load wild-type ecModel
load('ecModels/ecYeastGEM/model/ecYeastGEM_batch.mat')
ecModel = ecModel_batch;
clc
%% 3.- Set case-specific variables
%Indicate carbon source uptake reaction name
c_source   = 'D-glucose exchange (reversible)';
%Indicate carbon source mol. weight in grams/mmol
CS_MW      = 0.18015;
%Indicate the name of the growth (or biomass exchange) pseudoreaction
growth_rxn = 'growth';
%get current directory
current = pwd;
%product name provide a product name
product_name = '2-phenylethanol';
%results folder name
results_folder = '../../results/2_phenylethanol_targets';
mkdir(results_folder)
%% 4.- Generate production model
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
    disp(' ')
    if isempty(target_pos)
       error(['Exchange reaction for ' product_name ' could not be created'])
    else
        target_rxn = prod_ecModel.rxns(target_pos);
        prod_ecModel = setParam(prod_ecModel,'obj',target_pos,1);
        disp([ecModel.rxnNames{target_pos} ' has been set as cellular objective'])
    end
    disp(' ')
else
    error('Your provided product is not part of the used ecModel metabolites list')
end 

%% 5.-Set constraints
cd method
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
%Check if model can carry flux for the target rxn
flux = haveFlux(const_ecModel,1-12,target_pos);
if flux 
   disp(['* Your ecModel is suited for production of: ' product_name])
else
    error('Your ecModel is not suited for production, please revise your production pathway connectivity')
end
disp(' ')
%% 6.- Get max. biomass yield
WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
%Obtain a suboptimal yield value to run ecFactory
expYield = 0.49*WT_yield;
disp('* The ecFactory method will scan flux distributions spanning from')
disp(['a fixed biomass yield of: ' num2str(0.5*expYield) ' to: ' num2str(0.5*expYield) ' [g biomass/g carbon source]']);
%% 7.- Run ecFactory method
% Check compatibility with the method
const_ecModel = check_enzyme_fields(const_ecModel);
try
    [optStrain,candidates,step] = run_ecFactory(const_ecModel,target_rxn,const_ecModel.rxns(CS_index),expYield,CS_MW,results_folder);
   catch
    disp('The model is not suitable for robust ecFSEOF')
end
cd (current)
%%
