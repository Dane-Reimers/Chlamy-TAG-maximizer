% Evaluate & visualize the effects of nitrogen starvation on the metabolic
% pathways in iCre1355 under autotrophic conditions. Constraints are
% applied with experimentally derived RNA-seq data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346115/

% Init cobra, do not update
initCobraToolbox(false);

% Load in data
model = readCbModel('../data/models/iCre1355_mixo.xml');
load ../data/transcript_gene_map;
expression_table = readtable('expression_data.xlsx');

% Set the objective to TAG demand
model = changeObjective(model, {'DM_tag_c'});


% Map the expression data to the transcript names in iCre1355
tmp = expression_table.Au10_2(2:end-1);
transcript = cell(length(tmp),1);
for i=1:length(tmp)
    transcript{i} = transcript_gene_map(...
        strcmp(tmp(i,1),transcript_gene_map(:,2)),1);
end

% Create a table so that we can deal with duplicate transcripts without
% messing with paired arrays
expression_value = expression_table.Var6(2:end-1);
tmp_table = table(transcript, expression_value);
transcript_table = cell2table(cell(0,2), 'VariableNames', {'gene', 'value'});

% Iterate through every row of the table, count misses
% TODO almost half of the genes are misses...
misses = {};
for i=1:height(tmp_table)
    transcripts = tmp_table.transcript{i};
    % If row contains transcripts, port them over to the new table
    if ~isempty(transcripts)
        T = table;
        T.gene = transcripts;
        % Fill value col in table with the expression_value at i
        T.value(:) = tmp_table.expression_value(i);
        transcript_table = [transcript_table;T];
    else
        misses{end+1} = tmp{i};
    end
end
clear tmp;
clear tmp_table;
clear transcript;
clear expression_value;

% Convert to a struct
transcript_struct = table2struct(transcript_table, 'ToScalar', true);

% map expressions to reactions
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, transcript_struct);

% use expressionRxns to set hard bounds on the max flux through every
% reaction, don't set reactions with no expression data (E-Flux)
maxExpr = prctile(expressionRxns, 90);
for i=1:length(expressionRxns)
    if expressionRxns(i) ~= -1
        thisFlux = (expressionRxns(i)/maxExpr)*100;
        model.ub(i) = min(thisFlux, 100);
    else
        model.ub(i) = 100;
    end
    
end

% Use FVA to find the min and max flux through every reaction
[minFlux, maxFlux] = fluxVariability(model, 0.9);

% Find the max flux through our TAG demand reaction. This is always .1091,
% regardless of the transcriptomics data that we use.
maxTAGFlux = maxFlux(find(strcmp(model.rxns, 'DM_tag_c')))

threshold_lb = 0;
threshold_ub = 0.1;


% Experiment with different methods of extracting a context specific model.
% None of these yielded results like those in the paper that we got the
% model from.
gimme_model = GIMME(model, expressionRxns, 1000);
opt = optimizeCbModel(gimme_model);

tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub);
tissueModel2 = iMAT(model, expressionRxns, threshold_lb, threshold_ub, 1e-8, {'DM_tag_c'});
FBAsolution1 = optimizeCbModel(tissueModel);
FBAsolution2 = optimizeCbModel(tissueModel2);
[directionRxns, involvedMets, deadEnds] = draw_by_met(tissueModel2, {'tag16018111Z180[c]'}, 'true', 1, 'struc', {''}, FBAsolution1.x);
% tissueModel3 = iMAT(model, expressionRxns, threshold_lb, threshold_ub, 1e-8, {'TAG demand'});
% tissueModel4 = iMAT(model, expressionRxns, threshold_lb, threshold_ub, 1e-8, {'R_DM_tag_c'});
% tissueModel = changeObjective(tissueModel, 'DM_tag_c');


