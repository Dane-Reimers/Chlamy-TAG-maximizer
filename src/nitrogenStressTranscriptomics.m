% Evaluate & visualize the effects of nitrogen starvation on the metabolic
% pathways in iCre1355 under autotrophic conditions. Constraints are
% applied with experimentally derived RNA-seq data from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346115/

% Load in data
model = readCbModel('../data/models/iCre1355_auto.xml');
load ../data/transcript_gene_map;
expression_table = readtable('expression_data.xlsx');

% Set the objective to TAG demand and biomass
% TODO check if this is the problem? Maybe the model shouldn't have an
% objective before transcriptomics stuff??
model = changeObjective(model, {'DM_tag_c', 'Biomass_Chlamy_auto'});

% Map the expression data to the transcript names in iCre1355
% TODO this only reads in the 0 min col, we should read in all of them so
% we can compare pathways across the duration of the experiment
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

% Iterate through every row of the table
for i=1:height(tmp_table)
    transcripts = tmp_table.transcript{i};
    % If row contains transcripts, port them over to the new table
    if ~isempty(transcripts)
        T = table;
        T.gene = transcripts;
        % Fill value col in table with the expression_value at i
        T.value(:) = tmp_table.expression_value(i);
        transcript_table = [transcript_table;T];
    end
end
clear tmp_table;
clear transcript;
clear expression_value;

% Convert str array to char array
transcript_table.gene = char(transcript_table.gene);

% Convert to a struct
transcript_struct = table2struct(transcript_table, 'ToScalar', true);

% map expressions to reactions
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, transcript_struct);

threshold_lb = 0;
threshold_ub = 0.1;

tissueModel = iMAT(model, expressionRxns, threshold_lb, threshold_ub, 1e-8, {'DM_tag_c'});
tissueModel = changeObjective(tissueModel, 'DM_tag_c');


