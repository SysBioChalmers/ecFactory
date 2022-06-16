function ECCs = getECCs(model,target,perturbation,tolerance)
% getECCs
%   
% Compute enzyme control coefficients over a provided target reaction for
% all enzymes in an ecModel. 
%   
%   model           an ecModel MATLAB structure
%   target          rxn index for the target reaction
%   perturbation    relative magnitude for the perturbation on each enzyme
%                   activity (default = 1.01)
%   tolerance       (Optional) a numerical tolerance for enzyme activity
%                   changes (default = 1E-9)
%
%   Usage:  ECCs = getECCs(model,target,perturbation,tolerance)
% 
%   Ivan Domenzain. Last Edited 2020-12-17

if nargin<4
    tolerance = 1E-4;
    if nargin<3
        perturbation = 1.01;
    end
end
enzymes = model.enzymes;
%Relax bounds and set target reaction as objective to maximize
model.ub(target) = 1000;
model.lb(target) = 0;
model            = setParam(model,'obj',target,1);
base_sol         = solveLP(model,1);
current          = pwd;
ECCs             = zeros(length(enzymes),1);
maxKcats         = zeros(length(enzymes),1);
cd GECKO/geckomat/utilities
for i=1:length(enzymes)
    enzyme     = enzymes{i};
    temp_model = model;
    if ~isempty(enzyme)
        %Get enzyme information
        [~,rxnIdx] = getKcat(model,enzyme);
        enzPos     = find(strcmpi(model.mets,['prot_' enzyme]));
        enzUsRxn   = find(strcmpi(model.rxns,['draw_prot_' enzyme]));
        Kcats      = model.S(enzPos,rxnIdx);
        %Get basal enzyme usage (i-th enzyme) and relax its bounds
        E_usg = base_sol.x(enzUsRxn);
        temp_model.lb(enzUsRxn) = 0;
        temp_model.ub(enzUsRxn) = 1000;
        %Modify original Kcats for i-th enzyme
        temp_model.S(enzPos,rxnIdx) = Kcats./perturbation;      
        %Get mean activity for the i-th enzyme (considering all reaction in
        %which it might be present)
        A_i     = mean(Kcats.*E_usg);
        new_sol = solveLP(temp_model);
        if ~isempty(new_sol.f)
            %Modified enzyme activity (A_i = Kcat*perturbation*E_i)
            A_i_new  = mean(Kcats*perturbation.*new_sol.x(enzUsRxn));
            %Quantify change in the target reaction flux
            delta_V  = (new_sol.x(target) - base_sol.x(target))/base_sol.x(target);
            %Relative change in i-th enzyme activity
            delta_Ea = (A_i_new-A_i)/A_i;
            objCC = (1/(perturbation-1))*delta_V;
            if abs(delta_V) > tolerance%abs(delta_V)>=tolerance & abs(delta_Ea)>=tolerance   %& delta_Ea~=0
                %If solution was feasible then calculate the control coefficient
                %objCC       = delta_V./delta_Ea;%(perturbation-1);%delta_Ea;
                ECCs(i)     = objCC;
                maxKcats(i) = max(Kcats);
            end
        end
    end      
end
cd (current)
ECCs = table(enzymes,model.enzGenes,model.enzNames,ECCs,'VariableNames',{'enzymes' 'genes' 'names' 'CC'});
end