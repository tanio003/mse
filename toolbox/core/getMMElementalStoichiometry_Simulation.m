function [MMComposition] = getMMElementalStoichiometry_Simulation(model,BOF_coefs)
% computes the elemental composition of each MM component of biomass in
% terms of mmol gMM-1. Multiply by the BOF composition to get mmol gDW-1.
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


for k = 1:numel(crudeFractions)
    clear elementMatrix
    crudeFraction = crudeFractions{k};
    % get the MM value (g gDW-1)
    crudeFractionInd = find(strcmp(crudeFraction,model.mets));
    BOFind = find(strcmp('BIOMASSCRUDE',model.rxns));
    %crudeFractionOpt = abs(full(model.S(crudeFractionInd,BOFind)));
    crudeFractionOpt = BOF_coefs(k);

    % get the individual components (mmol gMM-1)
    MMRxn = find(strcmp(synthesisRxns{k},model.rxns));
    LHSMetInd = find(model.S(:,MMRxn)<0);
    LHSMets = model.mets(LHSMetInd);
    % get rid of water and atp
    excludeMets = [{'H2O'},{'ATP'}];
    [junk, junk, excludeMetsInd] = intersect(excludeMets,LHSMets);
    LHSMetInd(excludeMetsInd) = [];
    LHSMets(excludeMetsInd) = [];
    
    LHSCoefs = -full(model.S(LHSMetInd,MMRxn)); % mmol gMM-1
    
    % get elemental composistion of each component
    elementAbbrevs = [{'C'},{'H'},{'N'},{'O'},{'P'},{'S'}];
    
    
    for i = 1:numel(LHSMetInd)
        formula = model.metFormulas(LHSMetInd(i));
        [elements, useMat, exitFlag, MW]=parseFormulas(formula);
        % match elements
        [junk, junk2, elementInd] = intersect(elementAbbrevs,elements.abbrevs);
        elementMatrix(i,:) = useMat(elementInd);
    end
    
    % now compute total elemental composition for MM
    elementMM(k,:) = sum(elementMatrix .* repmat(LHSCoefs,1,size(elementMatrix,2)),1); % mmol gMM-1
    
    % and the elemental composition in biomass
    elementMM_DW(k,:) = elementMM(k,:).*crudeFractionOpt; % mmol gDW-1

end





% store output
MMComposition = struct;
MMComposition.MM = elementMM;
MMComposition.DW = elementMM_DW;
MMComposition.elementAbbrevs = elementAbbrevs;
MMComposition.crudeFractions = crudeFractions;


end
