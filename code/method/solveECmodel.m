function [mutSolution,flag] = solveECmodel(mutant,model,method,prots,tol)
% solveECmodel
%
% Function that receives a model and a mutant of it and gets a flux
% distribution according to the selected optimization method. 'pFBA'
% performs a first FBA optimization, then fixes the optimal value for the
% objective function and performs a minimization of the total sum of fluxes
% (minimization of total protein usage for ecModels) for getting rid of
% loops. MOMA performs a quadratic programming simulation, minimizing the
% euclidean distance from the Wild-type flux distribution.
% 
% mutant        mutant model 
% model         a wild-type model structure
% method        'MOMA' or 'pFBA'
% prots         Inxedes for the protein exchange reactions (ecModels)
% tol           (double) numerical tolerance for bounds fixation
%
% mutSolution   Optimal flux distribution for the mutant model
% flag          TRUE if sumulation was succesful
%
% Usage: [mutSolution,flag] = solveECmodel(mutant,model,method,prots)
%
% Last modified.  Ivan Domenzain 2020-07-16
%

if nargin<5
    tol = 1E-9;
end

minFlag = false;
if strcmpi(method,'MOMA')
    [mutSolution,~, flag] = qMOMA(mutant,model,1.01);
    disp(flag)
elseif strcmpi(method,'pFBA')
	%Get a first optimization
    mutSolution = solveLP(mutant);
    mutSolution = mutSolution.x;
    if ~isempty(mutSolution) && any(mutSolution)
        index = find(mutant.c);
        if ~isempty(prots)
            %Fix optimal value for the objective and then minimize the total
            %sum of protein usages
            mutant.lb(index) = (1-tol)*mutSolution(index);
            mutant.c(:)      = 0;
            mutant.c(prots)  = -1;
            mutSolution      = solveLP(mutant,1);
            if ~isempty(mutSolution.x) && any(mutSolution.x)
                mutSolution = mutSolution.x;
                minFlag    = true;
            end
        end
        %If protein usage minimization does not work or model is not EC, then 
        %apply simple pFBA
        if ~minFlag
            mutSolution = solveLP(mutant,1);
            mutSolution = mutSolution.x;
        end
    else 
        mutSolution = zeros(length(mutant.rxns),1);
    end
end

if ~isempty(mutSolution) && any(mutSolution)
    flag = 1;
else 
    flag = 0;
end

end