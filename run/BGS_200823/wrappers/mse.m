clear
close all

%% run MSE
% Runs MSE either in batch (SLURM array index) or locally
% This script serves as a wrapper for runMSE, and is currently set up to
% run on a cluster. 

% John R. Casey
% 20210221

%% Change directory
% change this path to where you have MSE installed
% rootPath = '/Users/jrcasey/Documents/MATLAB/GitHub/MSE_Standalone';
rootPath = '/Users/tatsurotanioka/Desktop/Project/mse';
runID = 'BGS_200823';                         
runPath = strcat(rootPath,'/run/',runID);

% make sure the toolbox is in Matlab's path
toolboxPath = '/Users/tatsurotanioka/Desktop/Project/mse/toolbox';
addpath(genpath(toolboxPath));
% Remove default data path
rmpath(strcat(rootPath,'/data'));

%% Load data
% These files are generated in BatchSetup.m
cd(runPath)
load('data/output/FileNames.mat')
load('data/output/Gridding.mat')
load('data/output/Options.mat')
load('data/output/CruiseData.mat')

%% Get SLURM index

% preallocate a 3d matrix of dimensions nStations, nZ, nStr
idxMat = zeros(Gridding.nStations,Gridding.nZ,Gridding.nStr);
nIterations = size(idxMat,1).*size(idxMat,2).*size(idxMat,3); % use this number for the SLURM array length

% get job array index (comment out if running the missing set; see below)
job_array_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));
job_array_idx = 7866;
% locate coordinates
[a,b,c] = ind2sub(size(idxMat),job_array_idx);
station_idx = a;
depth_idx = b;
strIdx = c;

%% Run missing files
% After running compile_AMT_subjobs.m, sometimes there are missing entries.
% Save those missing entry indices in a file called missingFileNo.mat and
% run those locally or back on the server from here: (comment out
% otherwise)
%  job_array_idx_temp = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%  load('data/output/missingFileNo.mat');
%  job_array_idx = missingFileNo(job_array_idx_temp);

%% Load strain model

strName = Gridding.strNameVec{strIdx};

load(strcat('data/GEM/StrMod/',strName,'.mat'));

%% Load enviromental data for indexed grid point

envDat = getEnvVars(CruiseData,station_idx,depth_idx);

%% To runMSE manually
model = mod; % 
%% Run MSE

tic
[Solution] = runMSE(mod,envDat,FileNames,Gridding,Options);
dt = toc;
Solution.runtime = dt;
Solution.strName = strName;
Solution.z = Gridding.depthVec(depth_idx);
Solution.station = Gridding.stationsVec{station_idx};

%% Save solution

% save local copy
%save(cell2str(FileNames.destination_fileName),'FullSolution');

% save on server
% save(strcat('/nobackup1/jrcasey/','Solution_',num2str(job_array_idx),'.mat'),'Solution');
save(strcat(runPath,'/data/output/Solution/Solution_',num2str(job_array_idx),'.mat'),'Solution');



