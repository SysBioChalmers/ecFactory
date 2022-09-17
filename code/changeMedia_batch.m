function [model,pos] = changeMedia_batch(model,c_source,media,flux)
%
%changeMedia_batch
% 
% Function that modifies the ecModel and makes it suitable for batch growth
% simulations on different carbon sources.
%
% INPUT:
%   - model:  An enzyme constrained model
%   - media:  Media type ('YEP' for complex, 'MAA' minimal with Aminoacids,
%             'Min' for minimal media)
%   - anox:   (optional) TRUE if anaerobic conditions are desired, DEFAULT=
%             FALSE
%   - flux:   (Optional) A cell array with measured uptake fluxes in mmol/gDwh
%
% OUTPUT:
%   - model: The ECmodel with the specified medium constraints
%
% Ivan Domenzain        2020-01-17


if nargin<3
    media = 'Min';
end


% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange (reversible)'

%first block any uptake
[rxnIDs,exchange]  = getExchangeRxns(model);
exchange           = exchange(find(contains(rxnIDs,'_REV')));
model.ub(exchange) = 0;
pos = getComponentIndexes(model,c_source);
%Block O2 and glucose production (avoids multiple solutions):
model.ub(strcmp(model.rxnNames,'oxygen exchange'))    = 0;
model.ub(strcmp(model.rxnNames,'D-glucose exchange')) = 0;
%Find substrate production rxn and block it:
pos_rev = strcmpi(model.rxnNames,c_source(1:strfind(c_source,' (reversible)')-1));
model.ub(pos_rev) = 0;
%For growth on fructose and mannose the transport takes place in a passive
%way. [Boles & Hollenberg, 2006]
if strcmp(c_source,'D-fructose exchange (reversible)')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1134')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1134')) = 0;
elseif strcmp(c_source,'D-mannose exchange (reversible)')
    model.S(strcmp(model.mets,'s_0796'),strcmp(model.rxns,'r_1139')) = 0;
    model.S(strcmp(model.mets,'s_0794'),strcmp(model.rxns,'r_1139')) = 0;
end
%The media will define which rxns to fix:
if strcmpi(media,'YEP')
    N = 25;     %Aminoacids + Nucleotides
elseif strcmpi(media,'MAA')
    N = 21;     %Aminoacids
elseif strcmpi(media,'Min')
    N = 1;      %Only the carbon source
end
%UB parameter (manually optimized for glucose on Min+AA):
b = 0.02;
%UB parameter (manually optimized for glucose complex media):
c = 2;
%Define fluxes in case of ec model:
if nargin < 4   %Limited protein    
    if N>1
       flux    = b*ones(1,N);
       if N>21
           flux(22:25) = c;
       end
    end
    flux(1) = 1;
end
%Fix values as UBs:
for i = 1:N
    model.ub(pos(i)) = flux(i);
end
%Allow uptake of essential components
model = setParam(model, 'ub', 'r_1654_REV', +1000); % 'ammonium exchange';
model = setParam(model, 'ub', 'r_2100_REV', +1000); % 'water exchange' ;
model = setParam(model, 'ub', 'r_1861_REV', +1000); % 'iron(2+) exchange';
model = setParam(model, 'ub', 'r_1992_REV', +1000); % 'oxygen exchange';
model = setParam(model, 'ub', 'r_2005_REV', +1000); % 'phosphate exchange';
model = setParam(model, 'ub', 'r_2060_REV', +1000); % 'sulphate exchange';
model = setParam(model, 'ub', 'r_1832_REV', +1000); % 'H+ exchange' ;
model = setParam(model, 'ub', 'r_4593_REV', +1000); % 'chloride exchange' ;
model = setParam(model, 'ub', 'r_4595_REV', +1000); % Mn(2+) exchange
model = setParam(model, 'ub', 'r_4596_REV', +1000); % Zn(2+ exchange
model = setParam(model, 'ub', 'r_4597_REV', +1000); % Mg(2+) exchange
model = setParam(model, 'ub', 'r_2049_REV', +1000); % sodium exchange
model = setParam(model, 'ub', 'r_4594_REV', +1000); % Cu(2+) exchange
model = setParam(model, 'ub', 'r_4600_REV', +1000); % Ca(2+) exchange
model = setParam(model, 'ub', 'r_2020_REV', +1000); % potassium exchange
%Block some production fluxes
model = setParam(model, 'ub', 'r_1663', 0); % bicarbonate exchange
model = setParam(model, 'ub', 'r_4062', 0); % lipid backbone exchange
model = setParam(model, 'ub', 'r_4064', 0); % lipid chain exchange
%Allow biomass production 
model = setParam(model, 'ub', 'r_2111', +1000); % growth
if ~strcmpi(media,'Min')
    model = setParam(model, 'ub', 'r_1654', 0); % 'ammonium exchange';
    model = setParam(model, 'ub', 'r_2091', 0); % 'ammonium exchange';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pos = getComponentIndexes(model,c_source)
    pos(1)  = find(strcmpi(model.rxnNames,c_source));
    pos(2)  = find(strcmpi(model.rxnNames,'L-alanine exchange (reversible)'));
    pos(3)  = find(strcmpi(model.rxnNames,'L-arginine exchange (reversible)'));
    pos(4)  = find(strcmpi(model.rxnNames,'L-asparagine exchange (reversible)'));
    pos(5)  = find(strcmpi(model.rxnNames,'L-aspartate exchange (reversible)'));
    pos(6)  = find(strcmpi(model.rxnNames,'L-cysteine exchange (reversible)'));
    pos(7)  = find(strcmpi(model.rxnNames,'L-glutamine exchange (reversible)'));
    pos(8)  = find(strcmpi(model.rxnNames,'L-glutamate exchange (reversible)'));
    pos(9)  = find(strcmpi(model.rxnNames,'L-glycine exchange (reversible)'));
    pos(10) = find(strcmpi(model.rxnNames,'L-histidine exchange (reversible)'));
    pos(11) = find(strcmpi(model.rxnNames,'L-isoleucine exchange (reversible)'));
    pos(12) = find(strcmpi(model.rxnNames,'L-leucine exchange (reversible)'));
    pos(13) = find(strcmpi(model.rxnNames,'L-lysine exchange (reversible)'));
    pos(14) = find(strcmpi(model.rxnNames,'L-methionine exchange (reversible)'));
    pos(15) = find(strcmpi(model.rxnNames,'L-phenylalanine exchange (reversible)'));
    pos(16) = find(strcmpi(model.rxnNames,'L-proline exchange (reversible)'));
    pos(17) = find(strcmpi(model.rxnNames,'L-serine exchange (reversible)'));
    pos(18) = find(strcmpi(model.rxnNames,'L-threonine exchange (reversible)'));
    pos(19) = find(strcmpi(model.rxnNames,'L-tryptophan exchange (reversible)'));
    pos(20) = find(strcmpi(model.rxnNames,'L-tyrosine exchange (reversible)'));
    pos(21) = find(strcmpi(model.rxnNames,'L-valine exchange (reversible)'));
    pos(22) = find(strcmpi(model.rxnNames,'2''-deoxyadenosine exchange (reversible)'));
    pos(23) = find(strcmpi(model.rxnNames,'2''-deoxyguanosine exchange (reversible)'));
    pos(24) = find(strcmpi(model.rxnNames,'thymidine exchange (reversible)'));
    pos(25) = find(strcmpi(model.rxnNames,'deoxycytidine exchange (reversible)'));
    pos(26) = find(strcmpi(model.rxnNames,'D-glucose exchange (reversible)'));
end
