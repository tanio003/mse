clear
close all

%% MSE Batch Setup
% Run this before running mse either locally or to HPC sending off to the server.

% Defines directories, filenames, gridding parameters, and various options
% to set up a batch run. Loads and formats environmental data to match the
% gridding.

% Steps required to set up a run:
% Filenames - define paths to data, models, databases, and outputs
% Gridding - define domain and resolution of the batch run
% Optimization options - define step sizes, tolerances, etc for the solvers
% Other options - there are a few other miscellaneous parameters to attend
% to

% Saves FileNames, Gridding, CruiseData, and Options structures to ~/data/output folder

% The AMT-13 cruise is included as an example, but of course in principle
% any dataset will do. 

% John R. Casey
% 20210221

% Tatsuro Tanioka 
% 20220823

% Saving Grid for Bio-GO-SHIP data at the surface

%% Basic Setup
% Modify these options for your particular run experiments

runID = 'BGS_200826c';         
% Sample specific GEM (y/n)?
% If yes, only the specific GEM created from the metagenome at a specific
% location will be used. If no, all the reference GEMs will be used instead
% of metagenome-derived GEMs
Sample_specific = {'y'};

% Ecotype abundance weighted parameterization (y/n)?
% If yes, the ecotype abundance weighted model paramters (Genome Size, Cell radius, Gc%, Pmax,
% and OGT) from observation will be used. If no, parameters based 
% for the most dominant ecotype will be used.
Ecotype_weighted = {'y'};

% Set the Cruises to include (keep those to include)
Cruises_select = {'I09N','P18','AE1319','AE1319','BVAL46','NH1418','AMT28','I07N','C13.5'};
% Cruises_select = {'I09N', 'I07N'}

% set the confidence model (90%, 95%, 99%)
% alpha_level = 'alpha90';
% alpha_level = 'alpha95';
alpha_level = 'alpha99';
%% Change directory
% change this path to where you have MSE installed
rootPath = '/Users/tatsurotanioka/Desktop/Project/mse';

% change this path to where data and wradppers for specific run is stored
runPath = strcat(rootPath,'/run/',runID);
cd(runPath)
addpath(genpath(runPath));
% make sure the toolbox is in Matlab's path
toolboxPath = '/Users/tatsurotanioka/Desktop/Project/mse/toolbox';
addpath(genpath(toolboxPath));

%% Version
version = strcat({'_v'},datestr(date,'yyyymmdd'))

%% Version control and paths (save a copy of this to the server)
% These files are required to run MSE

FileNames = struct;
% Current PanGEM
FileNames.PanGEM_Path = '/Users/tatsurotanioka/Desktop/Project/mse/data/GEM/PanGEM.mat';
% Organism database
FileNames.orgDB_Path = strcat(rootPath,'/data/db/orgDatabase.csv');
% Ecotype database
FileNames.ecotypeDB_Path = strcat(rootPath,'/data/db/ecotypeDatabase.csv');
% Current strain models
FileNames.StrMod_Path = strcat(rootPath,'/data/GEM/StrMod.mat');
% Current metagenome sample models
FileNames.SampleMod_Path = strcat(rootPath,'/data/BioGOSHIP/SampleMod_Merged_',alpha_level,'.mat');
% List of strains to be analyzed
FileNames.strainList_Path = strcat(rootPath,'/data/db/strainList.mat');
% List of metagenome samples to be analyzed
FileNames.sampleList_Path = strcat(rootPath,'/data/BioGOSHIP/sampleList.mat');
% Current OGTDat
FileNames.OGTDat_Path = strcat(rootPath,'/data/db/OGTDat.csv');
% Cruise data
FileNames.CruiseDB_filename = strcat(rootPath,'/data/BioGOSHIP/BioGOSHIP_Metadata_INPROGRESS_merged_mse.csv');
% HyperPro profiles
FileNames.IrrDat_fileName = strcat(rootPath,'/data/BioGOSHIP/IrrDat_BioGOSHIP.mat');
% TpDat path
FileNames.TpDat_fileName = strcat(rootPath,'/data/db/TpDat.csv');
% PhysOpt and PigOpt constraints path
FileNames.PhysOptPigOptConstraints_fileName = strcat(rootPath,'/data/db/BOFConstraints.csv');
% PigDB path
FileNames.PigDB_fileName = strcat(rootPath,'/data/db/AbsorptionDatabase.csv');

% Destination for solution and GEMs

mkdir(strcat(runPath,'/data'));
mkdir(strcat(runPath,'/data/GEM'));
mkdir(strcat(runPath,'/data/GEM/StrMod'));
mkdir(strcat(runPath,'/data/output'));
mkdir(strcat(runPath,'/data/output/Solution'));
%% Save and output options

Options = struct;
Options.saveSolution = false;
Options.saveStrMod = false;
if Sample_specific == "y"
    Options.samplespecific = true;
else
    Options.samplespecific = false;
end

if Ecotype_weighted == "y"
    Options.ecotypeweighted = true;
else
    Options.ecotypeweighted = false;
end

%% Optimization solver options

Options.maxIter_TpOpt = 1000; % maximum iterations for TpOpt algorithm
Options.maxIter_physOpt = 1000; % maximum iterations for PhysOpt algorithm
Options.stepTol_TpOpt = 1e-8; % step tolerance for TpOpt algorithm

%% Other options

Options.fmax = 0.085; % maximum coverage of cell surface area by transporters (proportion; see Casey and Follows, 2020 Plos CompBio)

%% Gridding (save a copy of this to the server)

Gridding = struct;
% Stations to include in batch run (select IN cruise)
metadata_GOSHIP = readtable(FileNames.CruiseDB_filename);
metadata_GOSHIP_select = metadata_GOSHIP(ismember(metadata_GOSHIP.Cruise, Cruises_select),:);
Gridding.stationsVec = metadata_GOSHIP_select.Station'; % station names
Gridding.stationsVec2 = find(ismember(metadata_GOSHIP.Cruise, Cruises_select)); % corresponding station indices (in CruiseData)
Gridding.nStations = numel(Gridding.stationsVec);
Gridding.Cruise = metadata_GOSHIP_select.Cruise;
% Depth (m)
% Gridding.minZ = 10; % minimum depth
% Gridding.maxZ = 200; % maximum depth
% Gridding.intervalZ = 10; % depth interval
% Gridding.depthVec = Gridding.minZ:Gridding.intervalZ:Gridding.maxZ;
% Gridding.nZ = numel(Gridding.depthVec);
Gridding.minZ = 7; % minimum depth
Gridding.maxZ = 7; % maximum depth
Gridding.intervalZ = 7; % depth interval
Gridding.depthVec = Gridding.minZ:Gridding.intervalZ:Gridding.maxZ;
Gridding.nZ = numel(Gridding.depthVec);
% Wavelength (nm)
Gridding.minLambda = 400; % minimum wavelength (nm)
Gridding.maxLambda = 686; % maximum wavelength (nm)
Gridding.bandwidth = 2; % bandwidth (nm)
Gridding.lambdaVec = Gridding.minLambda:Gridding.bandwidth:Gridding.maxLambda; % nm wavelength
%% Get metagenome GEMs to analyse
if Options.samplespecific
    load(FileNames.sampleList_Path);
    Gridding.strNameVec = Gridding.stationsVec';
else
    load(FileNames.strainList_Path);
    Gridding.strNameVec = strList;
end
Gridding.nStr = numel(Gridding.strNameVec);

%% Import environmental data

% Some notes:
% - nutrient concentrations should be in uM
% - temperature is in degrees C
% - Spectral downwelling irradiance should be in W m-2 (nm bandwidth)-1

%% Import Cruise Data
% Load and format nutrient concentrations, temperature, and other cruise
% variables

CruiseDB = readtable(FileNames.CruiseDB_filename,'Delimiter',',','ReadVariableNames',true);
CruiseData = getCruiseData(CruiseDB,Gridding.stationsVec, Gridding.depthVec);

%% Import spectral irradiance
% load and format HyperPro profiles
load(FileNames.IrrDat_fileName);
IrrDat = IrrDat_BioGOSHIP;
[IrrDat2] = standardizeIrr(IrrDat,Gridding.stationsVec, Gridding.lambdaVec,Gridding.depthVec); %mmoles photons m-2 h-1 bandwidth*nm-1
IrrDat3 = reshape([IrrDat2{:}],numel(Gridding.depthVec),numel(Gridding.lambdaVec),numel(IrrDat2));
CruiseData.IrrDat = IrrDat3;
CruiseData.PAR = squeeze(nansum(IrrDat3,2));

%% Parse sample GEM models and save to data/GEM/StrMod/ folder
load(FileNames.StrMod_Path);
load(FileNames.SampleMod_Path);
delete(strcat(runPath,'data/GEM/StrMod/*.mat'));
for a = 1:Gridding.nStr
     strName = Gridding.strNameVec{a};
     if Options.samplespecific
        if alpha_level == "alpha90"
            StrMod = SampleMod_Merged_alpha90;
        elseif alpha_level == "alpha95"
            StrMod = SampleMod_Merged_alpha95;
        elseif alpha_level == "alpha99"
            StrMod = SampleMod_Merged_alpha99;
        end
     end
     mod = StrMod.(strName);
     mod.id = strName;
     save(strcat(runPath,'/data/GEM/StrMod/',strName,'.mat'),'mod');
end

%% Save to output folder

save(strcat(runPath,'/data/output/Gridding.mat'),'Gridding');
save(strcat(runPath,'/data/output/FileNames.mat'),'FileNames');
save(strcat(runPath,'/data/output/Options.mat'),'Options');
save(strcat(runPath,'/data/output/CruiseData.mat'),'CruiseData');

