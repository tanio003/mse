function Solution = runMSE(model,envDat,FileNames,Gridding,Options)

% Run MSE pipeline for one strain at one grid point

% Inputs
% strName       -       String. StrainID (e.g., MED4)
% Gridding      -       Struct. Contains the ranges and resolutions for
%                       stations, depths, and wavelengths. 
% EnvDat        -       Struct. Environmental data (nutrients, temp, light)
% FileNames     -       Struct. Contains a list of filenames to be used for
%                       loading strain models, cruise data, and databases.
% Options       -       Struct. Contains various options like 'save',
%                       'maxIter', whether to save the acclimated StrMod's. 

% Outputs
% Solution      -       Struct. Solution contains growth rates, fluxes, BOF
%                       coefficients, uptake rate bounds, pigment specific
%                       absorption, transporter numbers, and cell size.


%% Prepare strain model for simulation
if Options.samplespecific
    StrMod2 = prepEnvMod_metagem(model, FileNames, Options);
    StrMod2.description = strcat('Un-acclimated metagenome-derieved model initialized with size and composition constraints - Prochlorococcus metagenome sample ID',{' '},model.id);
else
    StrMod2 = prepEnvMod(model,FileNames);
    StrMod2.description = strcat('Un-acclimated strain model initialized with size and composition constraints - Prochlorococcus strain',{' '},model.id);
end

%% Optimal Transport
% Load TpDat
TpDat = readtable(FileNames.TpDat_fileName,'ReadVariableNames',true,'Delimiter',',');

% Compile environmental data to retrieve and format TpOpt parameter table
TpDat = getTpOptParams(TpDat,envDat);

% Set up LP
[TpOptLP] = getTpOptLP(StrMod2, TpDat, Options);

% Solve LP
[n_opt,fval,exitflag,output,lambdaTemp,grad,hessian] = fmincon(TpOptLP);

% Calculate and store S_star
S_star = calc_S_star(StrMod2,n_opt(1),TpDat,envDat.T);

% Update model with TpRxns constraints
[StrMod3, vOut] = updateModel_TpOpt(StrMod2,n_opt,TpDat);

% Save optimal cell radius
r_opt = StrMod3.rInit;

%% MM Acclimation and Photoacclimation

% Load physOpt and pigOpt constraints together
Constraints = readtable(FileNames.PhysOptPigOptConstraints_fileName,'ReadVariableNames',true,'Delimiter',',');

% parse entries and get crude fractions
clb_ind = find(contains(Constraints.Constraint,'clb_'));
cub_ind = find(contains(Constraints.Constraint,'cub_'));
crudeFractions = strrep(Constraints.Constraint(clb_ind),'clb_','');

% Prepare pigOpt
% load PigDB
PigDB = readtable(FileNames.PigDB_fileName,'ReadVariableNames',true);

% Specify which pigments to include
PigsIncluded = [{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];
PigMW = [891.4731  905.4566  536.8726  568.8714];
nPigs = numel(PigsIncluded);

% Get reformated pigment specific absorption spectra
PigDat = getSpecificAbsorption(PigDB,PigsIncluded,PigMW,Gridding.lambdaVec); %m2 mmol-1 bandwidth(nm)^-^1

% Compute pigment absorption total
EdSpec = envDat.Ed;
a = nansum(PigDat .* repmat(squeeze(EdSpec)',1,nPigs) .* Gridding.bandwidth); % mmol photons [mmol pig]-1 h-1

% Set up LP
if Options.samplespecific
    [physOptLP] = getPhysOptLP_metagem(StrMod3, FileNames,Constraints, a, PigsIncluded, Options.maxIter_physOpt);
else
    [physOptLP] = getPhysOptLP(StrMod3,FileNames,Constraints, a, PigsIncluded, Options.maxIter_physOpt);
end
% Solve LP
[xOut,fval,exitflag,output] = fmincon(physOptLP);

% Save absorption fluxes
pigAbs = a.*xOut(12:15); % mmol photons gDW-1 h-1

% Update model
[StrMod4] = updateModel_pigPhysOpt(StrMod3,xOut,a,crudeFractions);
StrMod4.description = strcat('Acclimated strain model Prochlorococcus strain',{' '},model.id);

%% Temperature correction
if Options.samplespecific
    ecotypeDat = readtable(FileNames.ecotypeDB_Path,'Delimiter',',','ReadVariableNames',true);
    OGTDat = ecotypeDat(:, {'Ecotype', 'OGT'});
    T = [273.15 + envDat.T]; % convert to K
    TCorr = TemperatureCorr_metagem(OGTDat,model.id,T, FileNames, Options);
else
    OGTDat = readtable(FileNames.OGTDat_Path,'Delimiter',',','ReadVariableNames',true);
    T = [273.15 + envDat.T]; % convert to K
    TCorr = TemperatureCorr(OGTDat,model.id,T);
end
%% Solve FBA and temporarily store output
[sol hsSolOut] = solveLP(StrMod4,1);
if sol.stat==1
    fluxes = sol.x.*TCorr(1);
    growth = -sol.f.*TCorr(1);
    shadow = hsSolOut.y;
else
    fluxes = zeros(numel(StrMod4.rxns),1);
    growth = 0;
    shadow = zeros(numel(StrMod4.mets),1);
end
BOF = xOut;

%% Store and save output
% replace zero growth entries with nans
if growth < 1e-3
growth = NaN;
fluxes = nan(numel(fluxes),1);
BOF = nan(numel(BOF),1);
n_opt = nan(numel(n_opt),1);
r_opt = NaN;
pigAbs = nan(numel(pigAbs),1);
vOut = nan(numel(vOut),1);
S_star = nan(numel(S_star),1);
end

% Store in results structure
Solution = struct;
Solution.Growth = growth;
Solution.Fluxes = fluxes;
Solution.Shadow = shadow;
Solution.BOF_coefs = BOF;
Solution.TpOpt = n_opt;
Solution.r_opt = r_opt;
Solution.pigAbs = pigAbs;
Solution.uptakeBounds = vOut;
Solution.S_star = S_star;
if Options.saveStrMod
    Solution.StrMod = StrMod4;
end

% Store growth rates for each model version
sol1 = solveLP(model,1); % Unacclimated, unconditioned
if sol1.stat==1
    Solution.StrMod1_growth = -sol1.f;
else
    Solution.StrMod1_growth = NaN;
end
sol2 = solveLP(StrMod2,1); % Unacclimated, conditioned
if sol2.stat==1
    Solution.StrMod2_growth = -sol2.f;
else
    Solution.StrMod2_growth = NaN;
end
sol3 = solveLP(StrMod3,1); % Tp acclimated, conditioned
if sol1.stat==1
    Solution.StrMod3_growth = -sol3.f;
else
    Solution.StrMod3_growth = NaN;
end
sol4 = solveLP(StrMod4,1); % Tp and MM acclimated, conditioned
if sol1.stat==1
    Solution.StrMod4_growth = -sol4.f;
else
    Solution.StrMod4_growth = NaN;
end
sol5 = solveLP(StrMod4,1); % Tp and MM acclimated, Temperature corrected, conditioned
if sol5.stat==1
    Solution.StrMod5_growth = -sol5.f.*TCorr;
else
    Solution.StrMod5_growth = NaN;
end

end