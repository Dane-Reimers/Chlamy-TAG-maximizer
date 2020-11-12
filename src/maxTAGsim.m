% Find the flux through all of the TAG reactions

% initCobraToolbox(false) %don't update the toolbox
% changeCobraSolver ('gurobi', 'all', 1);
model = readCbModel('../data/models/iCre1355_auto.xml');
% [involvedMets, deadEnds] = draw_by_rxn(model, {'Biomass_Chlamy_auto'}, 'true');
% tags = model.mets(startsWith(model.mets,'tag'),:);
% tags;

% Get all of the reactions with TAG in the name
tag_rxns = model.rxns(startsWith(model.rxns, 'TAG'), :);
FBAsolution = optimizeCbModel(model);
FBAsolution.x;
% [Involved_mets, Dead_ends] = draw_by_rxn(model, tag_rxns, 'true', 'struc', {''}, {''}, FBAsolution.x);
% [directi onRxns, involvedMets, deadEnds] = draw_by_met(model, tags(1), 'true', 1, 'struc', {''}, FBAsolution.x);