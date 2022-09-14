% outer wrapper for solvePFBA

% load a model
fileName = 'data/GEM/pStrMod/pGEM_MED4.mat';
load(fileName);

% check that the model can grow
initSol = solveLP(model,1);

% compute optimal growth enzyme concentrations
kcat_conv = (1/6.022e23) .* 1000 .* 3600; % convert to mmol enzyme-1 h-1
E0_n = abs(initSol.x) ./ (model.Kcat .* kcat_conv); % n enzymes gDW-1 
E0 = E0_n .* (1/6.022e23) .* model.rxnEnzymeMW'; % g enzyme gDW-1

% proteome fraction of dry biomass
protFrac = nansum(E0); % 0.43!!! not bad!

% average molecular weight of the proteome
enz_molFrac = E0_n ./ nansum(E0_n); % mol enzyme mol proteome-1
proteomeMW = nansum(enz_molFrac .* model.rxnEnzymeMW'); % g mol-1

% now compute the amount of each KO. This will be the vector x to be
% optimized in OPFBA. 
for a = 1:numel(model.genes)
    % rxn indices for gene
    rxnIdx = find(model.rxnGeneMat(:,a));
    if ~isempty(rxnIdx)
    % enzyme [complex] concentration
    KO_n(a) = nansum(E0_n(rxnIdx)); % KO's gDW-1
    KO_g(a) = KO_n(a) .* (1/6.022e23) .* model.geneProductMW(a); % g KO gDW-1
    else
        KO_g(a) = 0;
    end
end

% check that the proteome fraction is the same computed this way
protFrac2 = nansum(KO_g);
tol = protFrac .* 0.01;
if abs(protFrac2 - protFrac) < tol
    disp('Complex and KO proteome fractions are equal to within 1%')
end


%% Setup PFBA problem

% assign to initial vector
x0 = KO_g;

% store x0 in model structure
model.x0 = x0;

% set bounds for x
lb = zeros(numel(model.genes),1);
ub_f = 3; % let's arbitrarily set the upper bound to be 10x the maximum flux value.
ub = repmat(ub_f,numel(model.genes),1); 

% set up decision matrix
A = ones(numel(model.genes),1);

% assign upper bound on proteome fraction
b = ub_f.*numel(model.genes); % maximum proteome g gDW-1

% set up optimization problem
options = optimoptions('fmincon','ConstraintTolerance',1e-3,'MaxIterations',1000);
prob = struct;
prob.x0 = ones(1,numel(model.genes)); % this time, we start with a vector of all ones
prob.objective = @(x)OPFBA_norm(x,model);
prob.Aineq = A';
prob.bineq = b;
prob.Aeq = [];
prob.beq = [];
prob.lb = lb;
prob.ub = ub;
prob.nonlcon = [];
prob.solver = 'fmincon';
prob.options = options;

[x,fval] = fmincon(prob);

% proteome fraction
protFrac_x = nansum(x.*model.x0);
% WOOOO finally!!!! 0.52 g gDW-1


%% Phenotype models

limMod = model;
limIdx = find(strcmp('IronEX',limMod.rxns))
% get optimal value
limOpt = initSol.x(limIdx);
% constrain bounds to half optimal
limMod.lb(limIdx) = 0.5.*limOpt;

% solve OPFBA
limProb = prob;
limProb.objective = @(x)OPFBA_norm(x,limMod);
[lim_x,lim_fval] = fmincon(limProb);

% get the amound of iron in each protein
ironIdx = find(strcmp('Fe2',model.mets));
protSynthIdx = find(strcmp('Enzyme Synthesis',model.subSystems));
nIron = full(model.S(ironIdx,protSynthIdx));
IronEnz_idx = find(model.S(ironIdx,protSynthIdx))


%% OLD method (without normalization)


% assign to initial vector
x0 = KO_g;

lb = zeros(numel(model.genes),1);
ub = repmat(0.3,numel(model.genes),1);

A = ones(numel(model.genes),1);
b = 0.5; % maximum proteome g gDW-1

% set up optimization problem
options = optimoptions('fmincon','ConstraintTolerance',1e-3,'MaxIterations',1000);
prob = struct;
prob.x0 = x0;
prob.objective = @(x)OPFBA(x,model);
prob.Aineq = A';
prob.bineq = b;
prob.Aeq = [];
prob.beq = [];
prob.lb = lb;
prob.ub = ub;
prob.nonlcon = [];
prob.solver = 'fmincon';
prob.options = options;

[x,fval] = fmincon(prob);


