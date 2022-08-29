function newMod = prepEnvMod_metagem(model, FileNames,Options)
%% Prepare strain model for environmental simulation
% Sets constraints and S matrix coefficients for use within a full
% simulation based on environmental data from either a cruise or from
% Darwin output. 


% These two steps will initialize a strain model to represent growth under
% optimal conditions.

% Inputs 
% model         -       Structure. Strain model. 
% orgDatabase   -       Table. Contains all necessary information to initialize strain model 

% Outputs
% newMod        -       Structure. Strain model ready for acclimation.

% John R. Casey
% 20191023

% Modified for metagenome-based model Tatsuro Tanioka
% 20220829

%% Prepare strain model for simulation
% Which strain
strName = model.id;

% Gather sample information
CruiseDatabase = readtable(FileNames.CruiseDB_filename,'ReadVariableNames',true,'Delimiter',',');
SampleInd = find(strcmp(CruiseDatabase.Station,strName));
SampleInd = SampleInd(1);

ecotype = CruiseDatabase.Ecotype{SampleInd};
ecotypeDatabase = readtable(FileNames.ecotypeDB_Path,'ReadVariableNames',true,'Delimiter',',');
ecotypeInd = find(strcmp(ecotypeDatabase.Ecotype,ecotype));

% Get ecotype abundance
if Options.ecotypeweighted
    HLI_frac = CruiseDatabase.HLI(SampleInd);
    HLII_frac = CruiseDatabase.HLII(SampleInd);
    HLIII_IV_frac = CruiseDatabase.HLIII_IV(SampleInd);
    LLI_frac = CruiseData.LLI(SampleInd);
    LLII_III_frac = CruiseData.LLII_III(SampleInd);
    LLIV_frac = CruiseData.LLIV(SampleInd);
    Abundace_table = [HLI_frac;HLII_frac;HLIII_IV_frac;LLI_frac;LLII_III_frac;LLIV_frac];
end

% Get cell size and weight
if Options.ecotypeweighted
    cell_r = 1e-6.*dot(Abundace_table,ecotypeDatabase.Cell_radius);  % m cell-1
else
    cell_r = 1e-6.*ecotypeDatabase.Cell_radius(ecotypeInd); % m cell-1
end
cellDW = (4/3).*pi().*((1e6.*cell_r).^3).*280e-15.*2; % gDW cell-1

% Get genome composition
if Options.ecotypeweighted
    GenomeSize = dot(Abundace_table, ecotypeDatabase.Size_Mb) .* 1e6 .* 2; % total number of bases (not pairs)
    GC = dot(Abundace_table,ecotypeDatabase.GCpct).* 1/100;
else
    GenomeSize = ecotypeDatabase.Size_Mb(ecotypeInd) .* 1e6 .* 2; % total number of bases (not pairs)
    GC = ecotypeDatabase.GCpct(ecotypeInd) .* 1/100;
end

nts = [{'dATP'},{'dGTP'},{'dCTP'},{'dTTP'}]; % dNTP's in DNA
nA = GenomeSize .* (1-GC) .* 0.5;
nG = GenomeSize .* GC .* 0.5;
nC = nG;
nT = nA;
avogadro = 6.022e23; % molecules mole-1
nNTs = [nA, nG, nC, nT]; % molecules cell-1
ntMW = [491.18,507.18,467.16,482.17] - 177.97; % g mol-1
ntMW_ave = sum(ntMW .* (nNTs ./ GenomeSize)); % average MW DNA base (g mol-1)
genomeMass = sum(nNTs .* (1/avogadro) .* ntMW); % g cell-1
genomeMassFrac = genomeMass ./ cellDW; % g gDW-1
gDNA = 1./ntMW_ave; % g DNA
coef_nt = 1000.* gDNA .* [0.5.*(1-GC) 0.5.*GC 0.5.*GC 0.5.*(1-GC)]; % mmol NTP g DNA-1

%% Step 1 - Set Pmax
% Set Pmax
if Options.ecotypeweighted
    cellPmax = dot(Abundace_table, ecotypeDatabase.Pmax); % fg cell-1 h-1
else
    cellPmax = ecotypeDatabase.Pmax(ecotypeInd); % fg cell-1 h-1
end
PmaxDW = 1000.* 1e-15 .* cellPmax .* (1/12.011) .* (1./cellDW); % 
model.ub(find(strcmp('R00024',model.rxns))) = PmaxDW;

%% Step 2 - Adjust NTP coefficients and DNA BOF coefficient
% Assign BOF coefficients
BOF_ind = find(contains(model.rxns,'BIOMASSCRUDE'));
DNA_ind = find(contains(model.mets,'DNA'));
%StrMod1.S(DNA_ind,BOF_ind) = -genomeMassFrac;
[model2, checkSum] = BOFadjust(model,'DNA',genomeMassFrac);
DNASynthesis_ind = find(contains(model.rxns,'DNASynthesis'));
for i = 1:numel(nts)
    [NTP_ind(i)] = find(strcmp(model2.mets,nts{i}));
end
PP_ind = find(contains(model2.mets,'Diphosphate'));
model2.S(NTP_ind,DNASynthesis_ind) = -coef_nt;
model2.S(PP_ind,DNASynthesis_ind) = sum(coef_nt);

%% Step 3 - Adjust GAM and NGAM
[GAM, NGAM, model3] = updateMaintenance(model2,cell_r,cellDW);

%% Add the initial cell size as a model structure field
if Options.ecotypeweighted
    model3.rInit = dot(Abundace_table, ecotypeDatabase.Cell_radius);
else
    model3.rInit = ecotypeDatabase.Cell_radius(ecotypeInd);
end
newMod = model3;

end



