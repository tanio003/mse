function dHc = getBiomassEnthalpy(model,BOF_coefs,EnthalpyDB)
%% Compute enthalpy of combustion of biomass given some MM composition

crudeFractions = [{'DNA'},{'Lipid'},{'Carbohydrate'},{'Protein'},...
    {'VitaCofactors'},{'RNA'},{'NB'},{'Ions'},{'BioPool'},...
    {'Pigments'},{'CellWall'}];
synthesisRxns = [{'DNASynthesis'},  
    {'LipidSynthesis'        },
    {'CarbohydrateSynthesis' },
    {'ProteinSynthesis'      },
    {'VitaCofactorsSynthesis'},
    {'RNASynthesis'          },
    {'NBSynthesis'           },
    {'IonSynthesis'          },
    {'BioPoolSynthesis'      },
    {'PigmentSynthesis'      },
    {'CellWallSynthesis'     }]';



for a = 1:numel(crudeFractions)
    crudeFraction = crudeFractions{a};
    % get the MM value (g gDW-1)
    crudeFractionInd = find(strcmp(crudeFraction,model.mets));
    BOFind = find(strcmp('BIOMASSCRUDE',model.rxns));
    %crudeFractionOpt = abs(full(model.S(crudeFractionInd,BOFind)));
    crudeFractionOpt = BOF_coefs(a);

    % get the individual components (mmol gMM-1)
    MMRxn = find(strcmp(synthesisRxns{a},model.rxns));
    LHSMetInd = find(model.S(:,MMRxn)<0);
    LHSMets = model.mets(LHSMetInd);
    LHSCoefs = -full(model.S(LHSMetInd,MMRxn)); % mmol gMM-1
    
    % get dHc for each LHSMets
    for b = 1:numel(LHSMets)
        idx = find(strcmp(LHSMets{b},EnthalpyDB.METS));
        LHSMets_dHc(b) = EnthalpyDB.dHc(idx); % KJ mol-1
    end
    dHc_MM(a) = nansum(LHSMets_dHc .* LHSCoefs' .* 1e-3); % KJ gMM-1
    clear LHSMets_dHc
end
% need to add on the pigments
pigs = [{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'},{'alpha_Carotene'},{'Zeaxanthin'}];
pigsMW = [891.47326, 905, 536.87264, 568.87144]; % g mol-1
for a = 1:numel(pigs)
    idx = find(strcmp(pigs{a},EnthalpyDB.METS));
    dHc_pigs(a) = EnthalpyDB.dHc(idx) .* (1./pigsMW(a)); % KJ gPig-1
end

dHc_MM2 = [dHc_MM dHc_pigs];

dHc = sum(dHc_MM2 .* BOF_coefs');

end

