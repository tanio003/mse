function [mu] = OPFBA(x,model)

% x is given in units of mmol gDW-1
x = reshape(x,numel(x),1);
%% Assign x vector coefficients to protein synthesis reaction
% get enzyme metabolite indices
enzyme_idx = find(contains(model.mets,'E_K'));

% get ATP, ADP, H2O, and Orthophosphate indices
costMets = [{'ATP'},{'ADP'},{'H2O'},{'Orthophosphate'}];
for a = 1:numel(costMets)
    costMets_idx(a) = find(strcmp(costMets{a},model.mets));
end

% get index of protein
prot_idx = find(strcmp('Protein',model.mets));

% get index of proteome synthesis reaction
ProtSynRxn_idx = find(strcmp('ProteomeSynthesis',model.rxns));

% get BOF index
BOFidx = find(strcmp('BIOMASSCRUDE',model.rxns));

% compute grams proteome aggregated from x
proteome_Mass = nansum(x .* model.geneProductMW .*(1./1000)); % g 


%% Assign bounds based on x and kcat

% convert Kcat units
k_cat = model.Kcat .* 3600; % mmol of reaction per mmol enzyme per hour

% tolerance for setting flux bound
tol = 0.1;

% threshold for setting flux bound
tau = 1e-6;

for a = 1:numel(x)
    % get reactions corresponding to KO
    rxnIdx = find(model.rxnGeneMat(:,a));
    nRxns = numel(rxnIdx);
    if nRxns > 0
        fluxBound = x(a) .* k_cat(rxnIdx); % mmol gDW-1 h-1
        for b = 1:nRxns
            if ~isnan(k_cat(rxnIdx(b)))
                if model.ub(rxnIdx(b)) > 0 & model.lb(rxnIdx(b)) == 0
                    model.ub(rxnIdx(b)) = max([fluxBound(b) .* (1+tol),tau]);
                elseif model.ub(rxnIdx(b)) == 0 & model.lb(rxnIdx(b)) < 0
                    model.lb(rxnIdx(b)) = min([-fluxBound(b) .* (1+tol),-tau]);
                elseif model.lb(rxnIdx(b)) < 0 & model.ub(rxnIdx(b)) > 0
                    model.lb(rxnIdx(b)) = min([-fluxBound(b) .* (1+tol),-tau]);
                    model.ub(rxnIdx(b)) = max([fluxBound(b) .* (1+tol),tau]);
                end
            end
        end
    end
end


%% Adjust proteome synthesis coefficients
% LHS protein coefficients
model.S(enzyme_idx,ProtSynRxn_idx) = -x;
% RHS proteome in grams
model.S(prot_idx,ProtSynRxn_idx) = proteome_Mass;

% ATP requirements for peptide polymerization
nAA_p = cellfun(@numel,model.geneProduct_sequence)'; % AA's per protein
nATP_p = 4.306 .* nAA_p; % ATP per protein, or mmol ATP per mmol protein
nATP_pCoef = x .* nATP_p; % mmol ATP per protein
nATP_total = nansum(nATP_pCoef); % mmol ATP per proteome

model.S(costMets_idx([1 3]),ProtSynRxn_idx) = -nATP_total; % ATP and water coefficients
model.S(costMets_idx([2 4]),ProtSynRxn_idx) = nATP_total; % ADP and phosphate coefficients

%% Adjust biomass reaction for new proteome fraction
[model, checkSum] = BOFadjust(model,'Protein',proteome_Mass);

% Re-assign boundary on RuBis Carboxylase flux
RuBisC_idx = find(strcmp('R00024',model.rxns));
model.ub(RuBisC_idx) = 4.7;

%% Solve model

sol = solveLP(model);
if sol.stat
    mu = sol.f;
else
    mu = 0;
end


end
