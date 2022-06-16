%Clone GECKO and go to the right branch
clear
current = pwd;
git('clone https://github.com/SysBioChalmers/GECKO.git')
cd GECKO
git checkout feat/add_FSEOF_utilities
git pull origin
cd ..
clc
success    = [];
fail       = [];
models     = [];
objectives = [];
steps      = [];
objIdx     = [];
noFlux     = [];
fileNames  = dir('../../ModelFiles/production_ecModels');
strain_conditions = readtable('../../data/strain_conditions.txt','Delimiter','\t');
mkdir('../../results/ECC')
for i=1:length(fileNames)
    cd (current)
    file = fileNames(i).name;
    if contains(file,'.mat')
        modelName = file(1:end-4);
        modelName = strrep(modelName,'_WBG','');
        disp([modelName ':  #' num2str(i)])
        if startsWith(modelName,'ec')
            modelName = lower(modelName(3:end));
            load(['../../ModelFiles/production_ecModels/' file])
            model       = check_enzyme_fields(model);
            targetIndex = find(model.c);
            targetRxn   = model.rxns(targetIndex);
            models      = [models;{file}];
            objectives  = [objectives;model.rxnNames(targetIndex)];
            objIdx      = [objIdx;targetIndex];
            disp(' ')
            %Set media conditions
            c_source  = 'D-glucose exchange (reversible)';
            tempModel = changeMedia_batch(model,c_source,'Min');
            CS_MW     = 0.18015;
            CS_index  = find(strcmpi(tempModel.rxnNames,c_source));
            growthPos = find(strcmpi(tempModel.rxnNames,'growth'));
            %Unconstrain growth
            tempModel = setParam(tempModel,'lb',growthPos,0);
            tempModel = setParam(tempModel,'ub',growthPos,1000);
            %Get biomass yield for a unit glucose uptake rate
            tempModel = setParam(tempModel,'obj',growthPos,1);
            tempModel = setParam(tempModel,'ub',CS_index,1);
            solution  = solveLP(tempModel,1);
            %Check if model can carry flux for the target rxn
            flux = haveFlux(tempModel,1-12,targetIndex);
            if ~flux | solution.x(growthPos)<=1e-4
                %find specific conditions and strain background
                idx = find(strcmpi(strain_conditions.Chemicals,modelName));
                if ~isempty(idx)
                    medium    = strain_conditions.Medium{idx};
                    c_source  = [strain_conditions.CarbonSource{idx} ' exchange (reversible)'];
                    tempModel = changeMedia_batch(tempModel,c_source,medium);
                    solution  = solveLP(tempModel,1);
                    %                     if ~isempty(strain_conditions.deletion(idx))
                    %                         genes2delete = strsplit(strain_conditions.deletion{idx},',');
                    %                         columns = cell(length(genes2delete),1);
                    %                         columns(:,1) = {0};
                    %                         modifications = [genes2delete' columns columns];
                    %                         %Search genes in model
                    %                         [iA,iB] = ismember(modifications(:,1),model.genes);
                    %                         modifications = modifications(iA,:);
                    %                         if ~isempty(modifications)
                    %                             %model = getMutant(model,modifications,[],true);
                    %                         end
                    %                     end
                end
            end
            WT_yield  = solution.x(growthPos)/(solution.x(CS_index)*CS_MW);
            disp(['The maximum biomass yield is ' num2str(WT_yield) '[g biomass/g carbon source]']);
            %reformat chemical name
            modelName = strrep(modelName,' ','_');
            modelName = strrep(modelName,'-','_');
            modelName = strrep(modelName,',','');
            modelName = strrep(modelName,'(','');
            modelName = strrep(modelName,')','');
            modelName = regexprep(modelName,'[0-9]','');
            modelName = strrep(modelName,'__','');
            modelName = lower(modelName);
            while startsWith(modelName,'_')
                modelName = modelName(2:end);
            end
            
            if flux & solution.x(growthPos)>0
                %Run Enzyme sensitivity analysis
                try
                    %Set suboptimal growth rate as lb and get protein pool
                    %usage for product maximization
                    modelECCs = setParam(tempModel,'lb',growthPos,0.5*solution.x(growthPos));
                    modelECCs = setParam(modelECCs,'obj',targetIndex,1);
                    solution  = solveLP(modelECCs,1);
                    PpoolIdx  = find(strcmpi(modelECCs.rxns,'prot_pool_exchange'));
                    Ppool     = solution.x(PpoolIdx);
                    %Unconstrain GUR and constrain protein pool with low
                    %usage value, obtain ECCs for low glucose levels
                    modelECCs = setParam(modelECCs,'ub',CS_index,1000);
                    modelECCs = setParam(modelECCs,'ub',PpoolIdx,Ppool);
                    %ECCs      = getECCs(modelECCs,targetIndex);
                    %restablish original protein pool UB and obtain ECCs
                    %for high glucose levels
                    modelECCs = setParam(modelECCs,'ub',PpoolIdx,tempModel.ub(PpoolIdx));
                    ECCs      = getECCs(modelECCs,targetIndex);
                    ECCs.Properties.VariableNames = {'enzymes' 'genes' 'names' 'FCC'};
                    writetable(ECCs,['../../results/ECC/' modelName '_ECCs.txt'],'QuoteStrings',false,'Delimiter','\t')
                catch
                   disp('The model is not suitable for ECC analysis')
                end 
                resultsFolder = ['../../results/production_targets/' modelName '_targets'];
                mkdir(resultsFolder)
                WT_yield = 0.48; %WT yield on glucose minimal media
                expYield = 0.49*WT_yield;
                try
                    [optStrain,candidates,step] = robust_ecFSEOF(tempModel,targetRxn,'r_1714_REV',expYield,CS_MW,resultsFolder);
                    success = [success; {modelName}];
                    candidates.FCC = zeros(height(candidates),1);
                    try
                        ECCs =readtable(['../../results/ECC/' modelName '_ECCs.txt'],'Delimiter','\t');
                        [iA,iB] = ismember(candidates.enzymes,ECCs.enzymes);
                        candidates.FCC(find(iA)) = ECCs.FCC(iB);
                        writetable(candidates,[resultsFolder '/compatible_genes_results.txt'],'Delimiter','\t','QuoteStrings',false);
                    catch
                        disp('No ECCs file available for this case')
                    end
                    %Generate transporters target file
                    rxnsTable     = readtable([resultsFolder '/rxnsResults_ecFSEOF.txt'],'Delimiter','\t');
                    transpTargets = getTransportTargets(rxnsTable,tempModel);
                    writetable(transpTargets,[resultsFolder '/transporter_targets.txt'],'Delimiter','\t','QuoteStrings',false);  
                catch
                    disp('The model is not suitable for robust ecFSEOF')
                    fail = [fail; {modelName}];
                end
                disp(' ')
            else
                noFlux = [noFlux; {modelName}];
                disp('The model is not able to carry flux through the target reaction with the imposed constraints')
            end
        end
    end
end
