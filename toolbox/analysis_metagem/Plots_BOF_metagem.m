%% Plot BOF compositions
% Contour plots of each BOF component for ecotypes and the population. 

%% Load data
rootPath = '/Users/tatsurotanioka/Desktop/Project/mse';
runID = 'BGS_221114a';                         
runPath = strcat(rootPath,'/run/',runID);
addpath(genpath(runPath));
load(strcat(runPath,'/data/output/FullSolution_L2.mat'))
FullSolution = FullSolution_L2;
load(strcat(runPath,'/data/output/Options.mat'))
if ~Options.samplespecific
    load(strcat(runPath,'/data/GEM/strainList.mat'));
end

%% Add mse core tool box to the path
toolboxPath = '/Users/tatsurotanioka/Desktop/Project/mse/toolbox/core';
addpath(genpath(toolboxPath));
%% Parse out Gridding, CruiseData, FileNames, and PanGEM from FullSolution
Gridding = FullSolution.Gridding;
CruiseData = FullSolution.CruiseData;
PanGEM = FullSolution.PanGEM;

%% Retrieve solutions for population at each sample point
[PopulationSolution] = parseBGSSolutions(FullSolution);
%% Gridded domain
x = CruiseData.Lat;
y = Gridding.depthVec;

%% Assign strains to ecotypes
% orgDatabase = readtable('GitHub/mse_AMT/data/db/orgDatabase.csv','Delimiter',',','ReadVariableNames',true);
% ecotypeList = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'}];
% ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];
% 
% for i = 1:Gridding.nStr
%     strainInd = find(strcmp(Gridding.strNameVec{i},orgDatabase.StrainName));
%     ecotype{i} = orgDatabase.Ecotype{strainInd};
%     ecotypeInd(i) = find(strcmp(ecotype{i},ecotypeList));
% end
% 
% for i = 1:numel(ecotypeList)
%     strEco_idx{i} = find(ecotypeInd == i);
% end
% ecotypeList2 = [{'HLI'},{'HLII'},{'LLI'},{'LLII/LLIII'},{'LLIV'}];

%% BOF components, compositions, and elemental quotas
% list components
BOF_components = [{'DNA'},{'Lipid'},{'Carbohydrate'},{'Protein'},{'VitaCofactors'},{'RNA'},{'NB'},{'Ions'},{'BioPool'},{'Pigments'},{'CellWall'},{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];

%% Plot a BOF component against a cruise variable
% 
mkdir(strcat(runPath,'/Figures/BOF'))
Cruise_names = Gridding.Cruise;
Nuc_acid = PopulationSolution.BOF_coefs(:,1,:) + PopulationSolution.BOF_coefs(:,6,:);
% --- Fig 1 = PAR vs Macros --- %
fig1 = figure;
subplot(2,2,1) % Plot Protein % against PAR
    x =  CruiseData.PAR'./3600*1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    ylabel('Protein Fraction')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Carb % against PAR
    x =  CruiseData.PAR'./3600*1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,3,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % ylim([0.2, 0.4])
    ylabel('Carbohydrate Fraction')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Lipid + Cell wall vs PAR    
    x =  CruiseData.PAR'./3600*1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
        squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    ylabel('Lipid + Cell wall Fraction')
    set(gca,'FontSize',12)
subplot(2,2,4) % Chla against PAR
    x =  CruiseData.PAR'./3600*1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,12,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('Chla Fraction')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)
fileName = strcat(runPath,'/Figures/BOF/BGS_PAR_Macrofrac.eps');
saveas(fig1,fileName,'epsc');

% --- Fig 2 = DIN vs Macros --- %
fig2 = figure;
tot_N = CruiseData.Ammonium' + ...
    CruiseData.NitratePlusNitrite';

subplot(2,2,1) % Plot Protein % against DIN
    x =  tot_N'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Protein % ')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Carb against DIN
    x =  tot_N'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,3,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Carb % ')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Lipid + Cell wall vs DIN
    x =  tot_N'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
       squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    % y = squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Lipid + Cell wall %')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,4) % Nucleic acid (RNA) vs DIN
    x =  tot_N'./1000;
    % y = squeeze(Nuc_acid)';
    y = squeeze(PopulationSolution.BOF_coefs(:,6,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('RNA %')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
fileName = strcat(runPath,'/Figures/BOF/BGS_DIN_Macrofrac.eps');   
saveas(fig2,fileName,'epsc');

% --- Fig 3 = DIP vs Macros --- %
fig3 = figure;

subplot(2,2,1) % Plot Protein % against DIP
    x =  CruiseData.Orthophosphate'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Protein % ')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot Carb % against DIP
    y = squeeze(PopulationSolution.BOF_coefs(:,3,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Carb % ')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)    
subplot(2,2,3) % Lipid + Cell wall vs DIP
    x =  CruiseData.Orthophosphate'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
        squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('Lipid + Cell wall %')
    xlabel('DIP (uM)')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Nucleic acid (RNA) vs DIP
    x =  CruiseData.Orthophosphate'./1000;
    y = squeeze(PopulationSolution.BOF_coefs(:,6,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'XScale','log')
    ylabel('RNA %')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)   
fileName = strcat(runPath,'/Figures/BOF/BGS_DIP_Macrofrac.eps');
saveas(fig3,fileName,'epsc');    

% --- Fig 4 = SST vs Macros --- 
fig4 = figure;
subplot(2,2,1) % Plot Protein % against Temp
    x =  CruiseData.T';
    y = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('Protein %')
    xlabel('Temperature')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot Carb % against Temp
    x =  CruiseData.T';
    y = squeeze(PopulationSolution.BOF_coefs(:,3,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('Carb %')
    xlabel('Temperature')
    set(gca,'FontSize',12)    
subplot(2,2,3) % Plot Lipid + Cell wall against Temp    
    x =  CruiseData.T';
    y = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
            squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('Lipid + Cell wall %')
    xlabel('Temperature')
    set(gca,'FontSize',12)  
subplot(2,2,4) % Plot Nucleic acid (RNA) vs Temp    
    x =  CruiseData.T';
    y = squeeze(PopulationSolution.BOF_coefs(:,6,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('RNA %')
    xlabel('Temperature')
    set(gca,'FontSize',12)  
fileName = strcat(runPath,'/Figures/BOF/BGS_Temp_Macrofrac.eps');
saveas(fig4,fileName,'epsc');      

%% Determine elemental quotas and enthalpy for each strain, ecotype, and population for
% each element. 

%%%%%%%%%%%%%%%% Takes about 45 minutes. %%%%%%%%%%%%%%%%%%

% load EnthalpyDB
% fileName = '/Users/jrcasey/Documents/MATLAB/CBIOMES/Data/Cheminformatics/Enthalpy_BOF/Enthalpy_BOF.csv'; 
% EnthalpyDB = readtable(fileName,'Delimiter',',','ReadVariableNames',true);

if ~Options.samplespecific
% all strains
     for j = 1:Gridding.nStations
         for k = 1:Gridding.nStr
             BOF_coefs = squeeze(PopulationSolution.BOF_coefs(j,1:11,k));
             [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
             Quota(j,:,k) = sum(MMComposition.DW,1); % mmol element gDW-1
             BOF_coefs2 = squeeze(PopulationSolution.BOF_coefs(j,:,k));
%             % Enthalpy(j,k) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ gDW-1
         end
     end
end

if Options.samplespecific
    % Metagenome-derived GEMs
    for k = 1:Gridding.nStations
        BOF_coefs = squeeze(PopulationSolution.BOF_coefs(1,1:11,k));
        [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
        QuotaPop(k,:) = sum(MMComposition.DW,1); % mmol element gDW-1
        BOF_coefs2 = squeeze(PopulationSolution.BOF_coefs(1,:,k));
        % Enthalpy(k) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ gDW-1
    end
end

% % ecotype winners
% for i = 1:Gridding.nZ
%     for j = 1:Gridding.nStations
%         for k = 1:numel(ecotypeList)
%             BOF_coefs = squeeze(EcotypeSolution.(ecotypeList{k}).BOF_coefs(i,j,1:11));
%             [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
%             QuotaEco(i,j,:,k) = sum(MMComposition.DW,1); % mmol element gDW-1
%             BOF_coefs2 = squeeze(EcotypeSolution.(ecotypeList{k}).BOF_coefs(i,j,:));
%             % EnthalpyEco(i,j,k) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ gDW-1
%         end
%     end
% end
% % population
% for i = 1:Gridding.nZ
%     for j = 1:Gridding.nStations
%             BOF_coefs = squeeze(PopulationSolution.BOF_coefs(i,j,1:11));
%             [MMComposition] = getMMElementalStoichiometry_Simulation(PanGEM,BOF_coefs);
%             QuotaPop(i,j,:) = sum(MMComposition.DW,1); % mmol element ml-1
%             BOF_coefs2 = squeeze(PopulationSolution.BOF_coefs(i,j,:));
%             % EnthalpyPop(i,j) = getBiomassEnthalpy(PanGEM,BOF_coefs2,EnthalpyDB); % KJ ml-1
%     end
% end


%% Calculating C:H:N:O:S:P (molar) and C:Chl (by mass) for population (metagenome-based) and for strain
if Options.samplespecific
    CN_Pop = squeeze(QuotaPop(:,1) ./ QuotaPop(:,3));
    CP_Pop = squeeze(QuotaPop(:,1) ./ QuotaPop(:,5));
    NP_Pop = squeeze(QuotaPop(:,3) ./ QuotaPop(:,5));
    OP_Pop = squeeze(QuotaPop(:,4) ./ QuotaPop(:,5));
    HP_Pop = squeeze(QuotaPop(:,2) ./ QuotaPop(:,5));
    SP_Pop = squeeze(QuotaPop(:,6) ./ QuotaPop(:,5));
    Cchl_Pop = (2.*1e-3.*12.011.* QuotaPop(:,1)) ./ (squeeze(PopulationSolution.BOF_coefs(1,12,:)) + squeeze(PopulationSolution.BOF_coefs(1,13,:))); 
else
    CN_Str = squeeze(Quota(:,1,:) ./ Quota(:,3,:));
    CP_Str = squeeze(Quota(:,1,:) ./ Quota(:,5,:));
    NP_Str = squeeze(Quota(:,3,:) ./ Quota(:,5,:));
    OP_Str = squeeze(Quota(:,4,:) ./ Quota(:,5,:));
    HP_Str = squeeze(Quota(:,2,:) ./ Quota(:,5,:));
    SP_Str = squeeze(Quota(:,6,:) ./ Quota(:,5,:));
    Cchl_Str = squeeze((2.*1e-3.*12.011.* Quota(:,1,:))) ./ (squeeze(PopulationSolution.BOF_coefs(:,12,:)) + squeeze(PopulationSolution.BOF_coefs(:,13,:))); 
end
%% Calculate the demands of aerobic remineralization (Paulmier et al., 2009, BG, 6(5) with Sulfur) for population (metagenome-based) and for strain
if Options.samplespecific
    O2P_Pop = (-1).* (CP_Pop + 0.25.*HP_Pop - 0.5.*OP_Pop - 0.75.*NP_Pop + 1.5.*SP_Pop + 1.25);
    r_O2C_Pop = (-1).* O2P_Pop.*(1.0./CP_Pop);
    rsum_O2C_Pop =  (-1).*(O2P_Pop - 2.0.*NP_Pop) .* (1.0./CP_Pop);
else
    O2P_Str = (-1).* (CP_Str + 0.25.*HP_Str - 0.5.*OP_Str - 0.75.*NP_Str + 1.5.*SP_Str + 1.25);
    r_O2C_Str = (-1).* O2P_Str.*(1.0./CP_Str);
    rsum_O2C_Str =  (-1).*(O2P_Str - 2.0.*NP_Str) .* (1.0./CP_Str);   
end
%% --- Plot C:N:P against env variables ---
% --- Fig 5 = PAR vs C:N:P:-O2 --- %
fig5 = figure;
subplot(2,2,1) % Plot C:P against PAR
    x =  CruiseData.PAR'./3600*1000;
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    ylabel('C:P')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot N:P against PAR
    x =  CruiseData.PAR'./3600*1000;
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    ylabel('N:P')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Plot C:N against PAR
    x =  CruiseData.PAR'./3600*1000;
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    ylabel('C:N')
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ against PAR
    x =  CruiseData.PAR'./3600*1000;
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    xlabel('PAR [\mu mol quanta m^-^2 s^-^1]')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_PAR.eps');
saveas(fig5,fileName,'epsc');    

% --- Fig 6 = DIN vs C:N:P:-O2 --- %
fig6 = figure;
subplot(2,2,1) % Plot C:P against DIN
    x =  tot_N';
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot N:P against DIN
    x =  tot_N';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Plot C:N against DIN
    x =  tot_N';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('DIN [\mu mol/kg]')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ against DIN
    x =  tot_N';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    xlabel('DIN [\mu mol/kg]')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_DIN.eps');
saveas(fig6,fileName,'epsc');  

% --- Fig 7 = DIP vs C:N:P:-O2 --- %
fig7 = figure;
subplot(2,2,1) % Plot C:P against DIP
    x = CruiseData.Orthophosphate';
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot N:P against DIP
    x = CruiseData.Orthophosphate';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Plot C:N against DIP
    x = CruiseData.Orthophosphate';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('DIP [\mu mol/kg]')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ against DIP
    x = CruiseData.Orthophosphate';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x./1000,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    set(gca,'XScale','log')
    xlabel('DIP [\mu mol/kg]')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_DIP.eps');
saveas(fig7,fileName,'epsc');  

% --- Fig 8 = SST vs C:N:P:-O2 --- %
fig8 = figure;
subplot(2,2,1) % Plot C:P against SST
    x = CruiseData.T';
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('Temperature [degC]')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot N:P against SST
    x = CruiseData.T';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('Temperature [degC]')
    set(gca,'FontSize',12)
subplot(2,2,3) % Plot C:N against SST
    x = CruiseData.T';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('Temperature [degC]')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ against SST
    x = CruiseData.T';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    xlabel('Temperature [degC]')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_Temp.eps');
saveas(fig8,fileName,'epsc');  

% --- Fig 9 = Latitude vs C:N:P:-O2 --- %
fig9 = figure;
subplot(2,2,1) % Plot C:P against Lat
    x = CruiseData.Lat';
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('Latitude')
    set(gca,'FontSize',12)
subplot(2,2,2) % Plot N:P against Lat
    x = CruiseData.Lat';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('Latitude')
    set(gca,'FontSize',12)
subplot(2,2,3) % Plot C:N against Lat
    x = CruiseData.Lat';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('Latitude')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ against Lat
    x = CruiseData.Lat';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    xlabel('Lat')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_Lat.eps');
saveas(fig9,fileName,'epsc');  

% --- Fig 10 = Macro vs C:N:P:-O2 --- %
fig10 = figure;
subplot(2,2,1) % C:P vs Nucleic acid (RNA)
    x = squeeze(PopulationSolution.BOF_coefs(:,6,:))'
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x(i,:).*100,y(i,:),'filled')
        hold on
    end
    ylabel('C:P')
    xlabel('% RNA')
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'FontSize',12)
subplot(2,2,2) % N:P vs Prot/Nucleic acid
    x = squeeze(PopulationSolution.BOF_coefs(:,4,:))'./squeeze(Nuc_acid)';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x(i,:),y(i,:),'filled')
        hold on
    end
    ylabel('N:P')
    xlabel('Protein:Nucleic acid')
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'FontSize',12)    
subplot(2,2,3) % C:N against % Protein
    x = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x(i,:).*100,y(i,:),'filled')
        hold on
    end
    ylabel('C:N')
    xlabel('% Protein')
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'FontSize',12)
subplot(2,2,4) %  RQ vs Lipid + Cell wall
    x = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
            squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end        
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled')
        hold on
    end
    ylabel('RQ')
    xlabel('Lipid + Cell wall %')
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    set(gca,'FontSize',12)    
fileName = strcat(runPath,'/Figures/BOF/Macro_vs_CNP.eps');
saveas(fig10,fileName,'epsc'); 

%% Relationship between cell size and macros and C:N:P
% --- Fig 11 = Cell size vs Macros --- %
fig11 = figure;
subplot(2,2,1) % Plot Protein % against cell size
    x = PopulationSolution.r_opt';
    y = squeeze(PopulationSolution.BOF_coefs(:,4,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    ylabel('Protein Fraction')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)
subplot(2,2,2) % Carb against cell size
    x = PopulationSolution.r_opt';
    y = squeeze(PopulationSolution.BOF_coefs(:,3,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % ylim([0.2, 0.4])
    ylabel('Carbohydrate Fraction')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)
subplot(2,2,3) % Lipid + Cell wall vs cell size  
    x = PopulationSolution.r_opt';
    y = squeeze(PopulationSolution.BOF_coefs(:,2,:))' + ...
        squeeze(PopulationSolution.BOF_coefs(:,11,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    xlabel('cell radius (\mu m)')
    ylabel('Lipid + Cell wall Fraction')
    set(gca,'FontSize',12)
subplot(2,2,4) % Chla against PAR
    x = PopulationSolution.r_opt';
    y = squeeze(PopulationSolution.BOF_coefs(:,12,:))';
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    ylabel('Chla Fraction')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)
fileName = strcat(runPath,'/Figures/BOF/BGS_Radius_Macrofrac.eps');
saveas(fig11,fileName,'epsc');

% --- Fig 12 = Cell size vs. CNP --- %
fig12 = figure;
subplot(2,2,1) % C:P
    x = PopulationSolution.r_opt';
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)
subplot(2,2,2) % N:P 
    x = PopulationSolution.r_opt';
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)
subplot(2,2,3) % C:N
    x = PopulationSolution.r_opt';
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('cell radius (\mu m)')
    set(gca,'FontSize',12)    
subplot(2,2,4) % Plot RQ 
    x = PopulationSolution.r_opt';
    if Options.samplespecific
        y = r_O2C_Pop';
    else
        y = r_O2C_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    xlabel('cell radius (\mu m)')
    ylabel('RQ')
    set(gca,'FontSize',12)      
fileName = strcat(runPath,'/Figures/BOF/CNP_vs_Radius.eps');
saveas(fig12,fileName,'epsc');  

%% Relationship between growth rate and macromolecules , C:N:P
% --- Fig 13 = Growth rate vs. CNP and macros--- %
fig13 = figure;
subplot(2,2,1) % C:P
    x = PopulationSolution.Growth;
    if Options.samplespecific
        y = CP_Pop';
    else
        y = CP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:P')
    xlabel('Growth rate [h^-^1]')
    set(gca,'FontSize',12)
subplot(2,2,2) % C:N
    x = PopulationSolution.Growth;
    if Options.samplespecific
        y = CN_Pop';
    else
        y = CN_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('C:N')
    xlabel('Growth rate [h^-^1]')
    set(gca,'FontSize',12)    
 subplot(2,2,3) % N:P
    x = PopulationSolution.Growth;
    if Options.samplespecific
        y = NP_Pop';
    else
        y = NP_Str';
    end
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('N:P')
    xlabel('Growth rate [h^-^1]')
    set(gca,'FontSize',12)      
    
subplot(2,2,4) % Nucleic acids (RNA)
    x = PopulationSolution.Growth;
    y = squeeze(PopulationSolution.BOF_coefs(:,6,:))'
    for i = 1:size(y,1);
        scatter(x,y(i,:),'filled');
        hold on 
    end
    if ~Options.samplespecific
        legend(strList,'Location','northeastoutside');
    end
    pointsize = 10;
    % set(gca,'XScale','log')
    ylabel('RNA %')
    xlabel('Growth rate [h^-^1]')
    set(gca,'FontSize',12)       

fileName = strcat(runPath,'/Figures/BOF/BGS_mu_CNP.eps');
saveas(fig13,fileName,'epsc');    

%% C:N vs DIN again, but with all strains, ecotypes, and population
% ecotypeColors = varycolor(numel(ecotypeList));
% fig = figure
% for i = 1:numel(ecotypeList)
%     for j = 1:numel(strEco_idx{i})
%         %x = CruiseData.Orthophosphate(Gridding.stationsVec2,:)';
%         x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
%         y = Quota(:,:,1,strEco_idx{i}(j)) ./ Quota(:,:,3,strEco_idx{i}(j));
%         plot(x,y,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
%         hold on
%     end
% end
% for i = 1:numel(ecotypeList)
%     x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
%     y = QuotaEco(:,:,1,i) ./ QuotaEco(:,:,3,i)
%     plot(x,y,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15)
%     hold on
% end
% x = CruiseData.Nitrate(Gridding.stationsVec2,:)' + CruiseData.Nitrite(Gridding.stationsVec2,:)' + CruiseData.Ammonia(Gridding.stationsVec2,:)';
% y = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
% plot(x,y,'.k','MarkerSize',15)
% %set(gca,'YScale','log')
% set(gca,'XScale','log')
% ylabel('C:N')
% %xlabel('E_d(474) / E_d(450)')
% xlabel('DIN [nM]')
% set(gca,'FontSize',20)
% 
% 
% saveas(fig,'/Users/jrcasey/Documents/New Structure/Projects/CBIOMES_Project/mse_AMT_Project/Figures/BOF/CN_vs_DIN_ALL.eps','epsc');



%% Depth Integrated C:N:P

% for a = 1:Gridding.nStations
%     C = QuotaPop(:,a,1) .* 1e6 .* 1e-3; % mol m-3
%     C_int(a) = trapz(Gridding.depthVec,C); % mol m-2
%     N = QuotaPop(:,a,3) .* 1e6 .* 1e-3; % mol m-3
%     N_int(a) = trapz(Gridding.depthVec,N); % mol m-2
%     P = QuotaPop(:,a,5) .* 1e6 .* 1e-3; % mol m-3
%     P_int(a) = trapz(Gridding.depthVec,P); % mol m-2
% end
% 
% cmap = [1 1 1; jet(256)];
% figure
% 
% % C:N
% subplot(3,1,1)
% imagesc(x,y,QuotaPop(:,:,1) ./ QuotaPop(:,:,3))
% hc = colorbar
% colormap(cmap)
% xlabel('Latitude')
% ylabel('Depth')
% ylabel(hc,'C:N [mol C (mol N)^-^1]')
% set(gca,'FontSize',20)
% 
% % C:P
% subplot(3,1,2)
% imagesc(x,y,QuotaPop(:,:,1) ./ QuotaPop(:,:,5))
% hc = colorbar
% colormap(cmap)
% xlabel('Latitude')
% ylabel('Depth')
% ylabel(hc,'C:P [mol C (mol P)^-^1]')
% set(gca,'FontSize',20)
% subplot(3,1,3)
% imagesc(x,y,QuotaPop(:,:,3) ./ QuotaPop(:,:,5))
% hc = colorbar
% colormap(cmap)
% xlabel('Latitude')
% ylabel('Depth')
% ylabel(hc,'N:P [mol N (mol P)^-^1]')
% set(gca,'FontSize',20)
% % C:Chl
% figure
% z = (2.*1e-3.*12.011.* QuotaPop(:,:,1)) ./ (PopulationSolution.BOF_coefs(:,:,12) + PopulationSolution.BOF_coefs(:,:,13)); 
% imagesc(x,y,z)
% hc = colorbar
% colormap(cmap)
% xlabel('Latitude')
% ylabel('Depth')
% ylabel(hc,'C:Chl [g C (g dv-Chl)^-^1]')
% set(gca,'FontSize',20)


%% Enthalpy
% ecotypes
% fig = figure
% for a = 1:numel(ecotypeList)
%     subplot(5,1,a)
%     imagesc(x,y,EnthalpyEco(:,:,a))
%     set(gca,'clim',[75 85])
%     hc = colorbar
%     ylabel(hc,strcat(ecotypeList2{a},' \DeltaH_c (KJ gDW^-^1)'))
%     colormap('jet')
%     xlabel('Latitude')
%     ylabel('Depth')
%     set(gca,'FontSize',20)
% end
% 
% Enthalpy2 = Enthalpy
% Enthalpy2(find(Enthalpy2==0)) = NaN
% nanmin(nanmin(nanmin(Enthalpy2)))
% nanmax(nanmax(nanmax(Enthalpy2)))


%% Carb vs lipid
% z = (EcotypeSolution.HLI.BOF_coefs(:,:,2) ) ./ EcotypeSolution.HLI.BOF_coefs(:,:,3);
% 
% figure
% for a = 1:numel(ecotypeList)
%     subplot(5,1,a)
%     z = (EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,2) ) ./ EcotypeSolution.(ecotypeList{a}).BOF_coefs(:,:,3);
%     imagesc(x,y,z)
%     colorbar
%     colormap('jet')
%     xlabel('Latitude')
%     ylabel('Depth')
%     set(gca,'FontSize',20)
%     
% end


%% C:N and C:P contours

% figure
% subplot(2,1,1)
% z = QuotaPop(:,:,1) ./ QuotaPop(:,:,3);
% 
% imagesc(x,y,z)
% hc=colorbar
% ylabel(hc,'C:N [mol C (mol N)^-^1]')
% colormap('jet')
% xlabel('Latitude')
% ylabel('Depth')
% set(gca,'FontSize',20)
% subplot(2,1,2)
% z = QuotaPop(:,:,1) ./ QuotaPop(:,:,5);
% imagesc(x,y,z)
% hc=colorbar
% ylabel(hc,'C:P [mol C (mol P)^-^1]')
% colormap('jet')
% xlabel('Latitude')
% ylabel('Depth')
% set(gca,'FontSize',20)



%% Cell Size

% figure
% z = (4/3) .* pi() .* PopulationSolution.r_opt.^3 .* 280; % fg C cell-1
% imagesc(x,y,z)
% hc=colorbar
% ylabel(hc,'Cell size (fg C cell^-^1')
% cmap = [1 1 1; jet(256)];
% colormap(cmap)
% xlabel('Latitude')
% ylabel('Depth')
% set(gca,'FontSize',20)


%% Chl b/a
% ecotypeColors = varycolor(numel(ecotypeList));
% load('/Users/jrcasey/Documents/MATLAB/GitHub/mse_AMT/data/envData/IrrDat.mat');
% aPeak_lambda = 450;
% bPeak_lambda = 474;
% aPeak_ind = find(IrrDat.Lambda == aPeak_lambda);
% bPeak_ind = find(IrrDat.Lambda == bPeak_lambda);
% for i = 1:Gridding.nStations
%     Ed_ratio(:,i) = IrrDat.Data{Gridding.stationsVec2(i)}(:,bPeak_ind) ./ IrrDat.Data{Gridding.stationsVec2(i)}(:,aPeak_ind);
% end
% 
% figure
% subplot(2,1,1)
% for i = 1:numel(ecotypeList)
%     for j = 1:numel(strEco_idx{i})
%         x1 =  Ed_ratio;
%         y1 = StrainSolution.BOF_coefs(:,:,13,strEco_idx{i}(j)) ./ StrainSolution.BOF_coefs(:,:,12,strEco_idx{i}(j));
%         plot(x1,y1,'.','MarkerFaceColor',ecotypeColors(i,:),'MarkerEdgeColor',ecotypeColors(i,:),'MarkerSize',5);
%         hold on
%     end
% end
% for i = 1:numel(ecotypeList)
%     x1 = Ed_ratio;
%     y1 = EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,13) ./ EcotypeSolution.(ecotypeList{i}).BOF_coefs(:,:,12);
%     h{i} = plot(x1,y1,'.','MarkerEdgeColor',ecotypeColors(i,:),'MarkerFaceColor',ecotypeColors(i,:),'MarkerSize',15);
%     hold on
% end
% x1 =  Ed_ratio;
% y1 = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12);
% h{6} = plot(x1,y1,'.k','MarkerSize',15)
% ylabel('Divinylchlorophyll b/a ratio')
% xlabel('E_d(474) / E_d(450)')
% xlim([1 3])
% set(gca,'FontSize',20)
% legend([h{1}(1) h{2}(1) h{3}(1) h{4}(1) h{5}(1) h{6}(a)],[ecotypeList2,{'Population'}])
% 
% subplot(2,1,2)
% z = PopulationSolution.BOF_coefs(:,:,13) ./ PopulationSolution.BOF_coefs(:,:,12)
% imagesc(x,y,z)
% hc=colorbar
% ylabel(hc,'Divinylchlorophyll b/a ratio')
% colormap('jet')
% xlabel('Latitude')
% ylabel('Depth')
% set(gca,'FontSize',20)






