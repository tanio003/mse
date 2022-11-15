% make_Irrdat_BioGOSHIP.m
% Script to make Irrdat matlab file for Bio-GO-SHP stations (surface data)
clear
%% Take a look at AMT_13 Irradiance matrix to see what variables are needed
load('~/Desktop/Project/JCasey_panGEM/mse/data/envData/IrrDat.mat');
IrrDat;
IrrDat.Lambda;  % wavelength array from 300 - 800 nm at 2 nm resolution
IrrDat.Depth;   % Depth array from 10 - 200 m at 10 resoultion (20 elements)
IrrDat.EdSurf;  % Incoming Spectral Irradiance at ocean surface for each wave band (251) at each station (54) in W m-2 nm-1
IrrDat.PAR;     % PAR at each station for each depth (20 x 54) in mol photons m-2 d-1
IrrDat.Station; % Station names (54)
IrrDat.Data;    % Spectral Irradiance Data at each depth (20) for each wave band (251) at each station(54) in W m-2 nm-1


%% Making Bio-GO-SHIP based IrrDat 
load('~/Desktop/Project/BioGOSHIP_MetaPanGEM/EnvData/OASIM_SpectralE_1deg_Monthly/BioGOSHIP_eds_ED_monthly_mse.mat');
load('~/Desktop/Project/BioGOSHIP_MetaPanGEM/EnvData/OASIM_PAR_1deg_Monthly/BioGOSHIP_PAR_OASIM_monthly_mse.mat');
metadata_GOSHIP = readtable('/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/EnvData/Metadata_merged/BioGOSHIP_Metadata_INPROGRESS_merged_mse.csv');
IrrDat_BioGOSHIP.Lambda = IrrDat.Lambda;
IrrDat_BioGOSHIP.Station = metadata_GOSHIP.Station';
IrrDat_BioGOSHIP.Depth  = [7];    % At 7 m only
IrrDat_BioGOSHIP.Data   = num2cell(eds_ED_Monthly_interp,2)'; % Spectral Downwelling irradiance in W m-2 nm-1 just below surface
IrrDat_BioGOSHIP.PAR    = PAR_Monthly';                       % PAR just below surface in mol photons m-2 d-1
save('IrrDat_BioGOSHIP.mat','IrrDat_BioGOSHIP');