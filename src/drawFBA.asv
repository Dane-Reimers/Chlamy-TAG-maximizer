% Load a model and draw it with Paint4Net

% initCobraToolbox(false) %don't update the toolbox
% changeCobraSolver ('gurobi', 'all', 1);
model = readCbModel('../data/models/iCre1355_auto.xml');
load ../data/transcript_gene_map;
% 'Biomass_Chlamy_auto'}
model = changeObjective(model, 'DM_tag_c');

% [involvedMets, deadEnds] = draw_by_rxn(model, {'Biomass_Chlamy_auto'}, 'true');
% tags = model.mets(startsWith(model.mets,'tag'),:);
% tags;

% Get all of the reactions with TAG in the name
% tag_rxns = model.rxns(startsWith(model.rxns, 'TAG'), :);
expression_table = readtable('expression_data.xlsx');

transcriptValue = expression_table.Var6(2:end-1);

expressionData = struct;
tmp = expression_table.Au10_2(2:end-1);
transcriptList = cell(length(tmp),1);
for i=1:length(tmp)
    transcriptList{i} = transcript_gene_map(...
        strcmp(tmp(i,1),transcript_gene_map(:,2)),1);
%     if transcriptList{i}
end
clear tmp;

newTransriptList = cell(length(transcriptList), 1);
newTransriptValue = zeros(length(transcriptValue), 1);
added = 0;
for i=1:length(transcriptList)
    
    if length(transcriptList{i}) > 1
        newTranscriptList{i + added} = char(transcriptList{i}{1});
        newTranscriptValue(i + added) = transcriptValue(i);
        for j=2:length(transcriptList{i})   
            added = added + 1;
            newTranscriptList{i + added} = char(transcriptList{i}{j});
            newTranscriptValue(i + added) = transcriptValue(i);
        end   
    else
        newTranscriptList{i + added} = char(transcriptList{i});
        newTranscriptValue(i + added) = transcriptValue(i);
    end
end

expressionData.gene = newTranscriptList;
expressionData.value = newTranscriptValue;
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData);
% Set the expression data for the demand reaction to be fully expressed so
% it doesn't get pruned by iMAT
expressionRxns(find(strcmp(model.rxns, 'DM_tag_c'))) = 1.0;
threshold_lb = 0;
threshold_ub = 0.1;

tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub, 1e-8, {});
tissueModel = changeObjective(tissueModel, 'DM_tag_c');
FBAsolution = optimizeCbModel(tissueModel);
no_stress = FBAsolution.f;

model = changeRxnBounds(model, 'EX_nh4_e', 0, 'l');
model = changeRxnBounds(model, 'EX_no3_e', 0, 'l');
FBAsolution = optimizeCbModel(model);
nitrogen_stress = FBAsolution.f;
[Involved_mets, Dead_ends] = draw_by_rxn(model, {'DM_tag_c'}, 'true', 'struc', {''}, {''}, FBAsolution.x);
% [directionRxns, involvedMets, deadEnds] = draw_by_met(model, {'tag16018111Z180[c]'}, 'true', 1, 'struc', {''}, FBAsolution.x);