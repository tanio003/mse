%% Post processing on local machine for Bio-GO-SHIP cruises

%% Directory

rootPath = '/Users/tatsurotanioka/Desktop/Project/mse';
runID = 'BGS_221114a';                         
runPath = strcat(rootPath,'/run/',runID);
addpath(genpath(runPath));

% make sure the toolbox is in Matlab's path
toolboxPath = '/Users/tatsurotanioka/Desktop/Project/mse/toolbox';
addpath(genpath(toolboxPath));

% Remove default data and wrapper path
rmpath(strcat(rootPath,'/data'));
rmpath(strcat(toolboxPath,'/wrappers'));

%% Load data
% These files are generated in BatchSetup.m
cd(runPath)
load('data/output/FileNames.mat')
load('data/output/Gridding.mat')
load('data/output/Options.mat')
load('data/output/CruiseData.mat')
load(strcat(rootPath,'/data/GEM/PanGEM.mat'))
%% Point to solutions directory
% ResultsDirectory = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Environmental_Data/Cruises/AMT13/Solution_20200914/';
ResultsDirectory = strcat(runPath,'/data/output/Solution/');

%% Compile results (takes less than a minute)
[FullSolution_L1] = get_BGS_Results_server(ResultsDirectory,Options,PanGEM,FileNames,Gridding,CruiseData);

% Save compiled results locally. 
%%%%%%%%%%%%%%%%%%% 
% IMPORTANT: DO NOT COMMIT FullSolution_L1.mat or FullSolution_L2.mat or you will
% overrun the github file size limit.
%%%%%%%%%%%%%%%%%%% 
%%
save(strcat(runPath,'/data/output/FullSolution_L1.mat'),'FullSolution_L1');