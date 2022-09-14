function [mu] = solvePFBA(x,model)

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

% compute average molecular weight of proteome
proteome_MW = nanmean(x .* model.geneProductMW); % g mol-1

% compute the proteome fraction of dry weight
protFrac_g = nansum(x); % g gDW-1

% compute moles of proteome per gDW
protFrac_mol = protFrac_g .* (1./proteome_MW); % mol gDW-1

% compute polymerization ATP costs
partFrac = x./protFrac_g; % relative proportion of each KO
nAA = cellfun(@numel,model.geneProduct_sequence); % AA's per KO
nAA_prot = nanmean(nAA); % average AA's per mole protein
nATP = 4.306 .* nAA_prot; % moles ATP per mole protein
nATP_gProt = nATP .* (1/proteome_MW) .* protFrac_g;

% assign x, proteome fraction, and costs to proteome synthesis reaction
model.S(enzyme_idx,ProtSynRxn_idx) = -x; % enzyme coefficients
model.S(costMets_idx([1 3]),ProtSynRxn_idx) = -nATP_gProt; % ATP and water coefficients
model.S(costMets_idx([2 4]),ProtSynRxn_idx) = nATP_gProt; % ADP and phosphate coefficients
model.S(prot_idx,ProtSynRxn_idx) = protFrac_g;
model.ub(ProtSynRxn_idx) = 1000;

% get index of MM protein synthesis reaction, we will need to block this
% reaction
MMProtSynRxn_idx = find(strcmp('ProteinSynthesis',model.rxns));
model.ub(MMProtSynRxn_idx) = 0;

%% Assign bounds based on x and kcat

% convert Kcat units (molecules enzyme-1 s-1) to FBA units (mmol g enzyme-1
% h-1)
k_cat = model.Kcat .* (1/6.022e23) .* 1000 .* 3600; % mmol g enzyme-1

% tolerance for setting flux bound
tol = 0.1;

for a = 1:numel(x)
    % get reactions corresponding to KO
    rxnIdx = find(model.rxnGeneMat(:,a));
    nRxns = numel(rxnIdx);
    if nRxns > 0
        fluxBound = x(a) .* k_cat(rxnIdx) .* 6.022e23 .* (1./model.geneProductMW(a)); % mmol gDW-1 h-1
        for b = 1:nRxns
            if model.ub(rxnIdx(b)) > 0 & model.lb(rxnIdx(b)) == 0
                model.ub(rxnIdx(b)) = fluxBound(b) .* (1+tol);
            elseif model.ub(rxnIdx(b)) == 0 & model.lb(rxnIdx(b)) < 0
                model.lb(rxnIdx(b)) = -fluxBound(b) .* (1+tol);
            elseif model.lb(rxnIdx(b)) < 0 & model.ub(rxnIdx(b)) > 0
                model.lb(rxnIdx(b)) = -fluxBound(b) .* (1+tol);
                model.ub(rxnIdx(b)) = fluxBound(b) .* (1+tol);
            end
        end
    end
end

%% Adjust biomass reaction for new proteome fraction
[model, checkSum] = BOFadjust(model,'Protein',protFrac_g);

%% Solve model

sol = solveLP(model);
if sol.stat
    mu = sol.f;
else
    mu = 0;
end


end
