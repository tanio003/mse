function TCorr = TemperatureCorr_metagem(OGTDat,sampleName,T,FileNames,Options)
% compute the temperature correction based on OGT and the in situ
% temperature

% Inputs
% OGTDat        -       Table. Optimal growth temperature database (based
%                       on the machine learning algorithm)
% StrainName    -       String. Strain ID
% T             -       Double. In situ growth temperature (K)

% Outputs
% TCorr         -       Double. Temperature correction. This value
% currently scales with the activation energy calculated for
% Prochlorococcus. No consideration of the temperature range of growth is
% made currently. Instead, the growth correction is forced to 1.0 for all
% in situ temperatures above OGT;

% John R. Casey
% 20190923

% Modified for metagenome-based models by Tatsuro Tanioka
% 20220829

%% Sample OGT
% Get ecotype abundance and calculate weighted mean OGT
CruiseDatabase = readtable(FileNames.CruiseDB_filename,'ReadVariableNames',true,'Delimiter',',');
SampleInd = find(strcmp(CruiseDatabase.Station,sampleName));
SampleInd = SampleInd(1);
ecotype = CruiseDatabase.Ecotype{SampleInd};
ecotypeDatabase = readtable(FileNames.ecotypeDB_Path,'ReadVariableNames',true,'Delimiter',',');
ecotypeInd = find(strcmp(ecotypeDatabase.Ecotype,ecotype));

if Options.ecotypeweighted
    HLI_frac = CruiseDatabase.HLI(SampleInd);
    HLII_frac = CruiseDatabase.HLII(SampleInd);
    HLIII_IV_frac = CruiseDatabase.HLIII_IV(SampleInd);
    LLI_frac = CruiseDatabase.LLI(SampleInd);
    LLII_III_frac = CruiseDatabase.LLII_III(SampleInd);
    LLIV_frac = CruiseDatabase.LLIV(SampleInd);
    Abundace_table = [HLI_frac;HLII_frac;HLIII_IV_frac;LLI_frac;LLII_III_frac;LLIV_frac];
    OGT = 273.15 + dot(Abundace_table, OGTDat.OGT);  % K
else
    OGT = 273.15 + OGTDat.OGT(SampleInd); % K
end
%% Arrhenius parameters

Ea = 5.2733e4; %J mol-1
R = 8.314;
A = 1;
OGT_rate = exp(-Ea./(R.*(OGT)));
%InSitu_rate = exp(-Ea./(R.*(T)));
InSitu_rate = exp(-Ea./(R.*(T))).*(1-exp(T-(OGT+2)));

% Lower bound linear correction (optional, comment out for just the
% Arrhenius function). Currently using 10 degrees less than OGT as the
% upper bound of this linear function, and 15 degrees less  than OGT as the
% lower bound (zero rate). 

if T < OGT - 10
    OGT_minus10_rate = exp(-Ea./(R.*(OGT-10)));
    LB_slope = OGT_minus10_rate / 5;
    dT = OGT - 10 - T;
    InSitu_rate = InSitu_rate - dT.*LB_slope;
end

% Check that we have a positive rate
if InSitu_rate < 0
    InSitu_rate = 0;
end


% temperature correction
TCorr = InSitu_rate ./ OGT_rate;


end