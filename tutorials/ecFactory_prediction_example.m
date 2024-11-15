% Prediction of metabolic engineering targets for increased production of 2-phenyethanol in Saccharomyces cerevisiae using ecYeastGEM and ecFactory

% 1.Get and load ecYeastGEM (Note: this model is an ecModel for S. cerevisiae. If the user wants to predict targets for other host organisms, they can load other ecModels as needed.)
git clone --quiet --depth=1 https://github.com/SysBioChalmers/ecModels.git
clc
load('ecModels/ecYeastGEM/model/ecYeastGEM_batch.mat')
ecModel = ecModel_batch;

% (Optional) If the product is heterologous product, you need to build product-specific ecModel by add the heterologous biosynthesis pathway of your product.
% take heme as an example
% ecModel = getHeme_ecYeastGEM(ecModel);


% 2. Clone GECKO inside ecFactory/code. 
% Note!! The ecFactory is currently only compatible with GECKO 2.0. 
% Adjustments made to the ecmodel structure in GECKO 3.0 are temporarily incompatible with ecFactory.
cd ../code
git clone --quiet --branch v2.0.3 --depth=1 https://github.com/SysBioChalmers/GECKO.git
addpath(genpath('GECKO'));


% 3.- Set case-specific parameters
% Create a structure with the following model-specific parameters:
% modelParam.CSrxn = preferred carbon source uptake reaction ID
% modelParam.CS_MW = Molecular weight of the preferred carbon source
% modelParam.growthRxn = Growth reaction name
% modelParam.rxnTarget = ID for the production target reaction

%Find carbon source uptake reaction 
CSname = 'D-glucose exchange (reversible)';
modelParam.CSrxn  = ecModel.rxns{strcmpi(ecModel.rxnNames,CSname)};
%Indicate carbon source mol. weight in grams/mmol
modelParam.CS_MW  = 0.18015;
%Indicate the name of the growth or biomass pseudoreaction
modelParam.growthRxn = ecModel.rxns{strcmpi(ecModel.rxnNames,'biomass pseudoreaction')};
%identify the production target reaction
product_name         = '2-phenylethanol';
target_pos           = find(strcmpi(ecModel.rxnNames,[product_name ' exchange']));
modelParam.rxnTarget = ecModel.rxns(target_pos);
disp(['*The production reaction for ' product_name ' has been found under the ID: ' modelParam.rxnTarget{1}])

% Provide a directory path for the results files
%results folder name
results_folder = '../tutorials/results/2_phenylethanol_targets';
mkdir(results_folder)
%get current directory
current = pwd;


% 4. Constrain ecModel
%Set media conditions for batch growth
const_ecModel = changeMedia_batch(ecModel,CSname,'Min'); %model-specific script

CS_index      = find(strcmpi(const_ecModel.rxns,modelParam.CSrxn));
growthPos     = find(strcmpi(const_ecModel.rxns,modelParam.growthRxn));
%Enable cellular growth
const_ecModel = setParam(const_ecModel,'lb',growthPos,0);
const_ecModel = setParam(const_ecModel,'ub',growthPos,1000);
%block artificial growth
revIndex = find(strcmpi(const_ecModel.rxnNames,'growth (reversible)'));
const_ecModel.lb(revIndex) = 0;
const_ecModel.ub(revIndex) = 0;
%Set a fixed unit glucose uptake rate
const_ecModel = setParam(const_ecModel,'ub',CS_index,1);


% 5. Get a feasible biomass yield range to be scanned by the ecFactory method
%Get biomass yield for a unit glucose uptake rate
const_ecModel = setParam(const_ecModel,'obj',growthPos,1);
solution      = solveLP(const_ecModel,1);
WT_yield      = solution.x(growthPos)/(solution.x(CS_index)*modelParam.CS_MW);
disp(['* The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);

%Obtain a suboptimal yield value to run ecFactory
expYield = 0.49*WT_yield;
disp('* The ecFactory method will scan flux distributions spanning from')

disp(['a suboptimal biomass yield of: ' num2str(0.5*expYield) ' to: ' num2str(2*expYield) ' [g biomass/g carbon source]']);


% 7. Run ecFactory method
try
    [optStrain,candidates,step] = run_ecFactory(const_ecModel,modelParam,expYield,results_folder,true);
catch
    disp('The model is not suitable for the ecFactory method')
end


% 8. Get a flux distribution of the optimal strain model
solution = solveLP(optStrain,1);
printFluxes(optStrain,solution.x,true)

