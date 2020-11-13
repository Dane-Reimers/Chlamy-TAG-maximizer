%
% This script performs a single gene analysis of the given models via 
% singleGeneDeletion method of the COBRA toolbox
%
% PREREQUISITES:
%   - have COBRA toolbox installed and initialized (run initCobraToolbox)
%
% Output description:
%
% each singleKOsim.{auto|mixo|hetero} contains the following data (see 
% also singleGeneDeletion.m from the COBRA toolbox)
%
% grRatio       Computed growth rate ratio between deletion strain and wild type
% grRateKO      Deletion strain growth rates (1/h)
% grRateWT      Wild type growth rate (1/h)
% hasEffect     Does a gene deletion affect anything (i.e. are any reactions
%               removed from the model)
% delRxns       List of deleted reactions for each gene KO
% fluxSolution  FBA/MOMA/lMOMA fluxes for KO strains
changeCobraSolver();
model_mixo = readCbModel('../data/models/iCre1355_mixo.xml');
load ../data/transcript_gene_map

model_mixo = changeObjective(model_mixo, "Biomass_Chlamy_mixo", 1);

printConstraints(model_mixo, -1000, 1000)

% process list of KO candidate genes via transcript_gene_map
tmp = unique(transcript_gene_map(2:end,2));
geneKOlist = cell(length(tmp),1);
for i=1:length(tmp)
    geneKOlist{i} = transcript_gene_map(...
        strcmp(tmp(i,1),transcript_gene_map(:,2)),1);
end
clear tmp

newGeneKOlist = cell(length(geneKOlist), 1);
added = 0;
for i=1:length(geneKOlist)
    
    if length(geneKOlist{i}) > 1
        newGeneKOlist{i + added} = char(geneKOlist{i}{1});
        for j=2:length(geneKOlist{i})   
            added = added + 1;
            newGeneKOlist{i + added} = char(geneKOlist{i}{j});
        end   
    else
        newGeneKOlist{i + added} = char(geneKOlist{i});
    end
end

[grRatioDble, grRateKO, grRateWT] = doubleGeneDeletion(model_mixo, 'FBA', newGeneKOlist(1:100));
for i = 1:length(grRatioDble)
    if grRatioDble(i) > 1.000001
        disp(grRatioDble(i));
    end
end

[singleKOsim.mixo.grRatio,...
    singleKOsim.mixo.grRateKO,...
    singleKOsim.mixo.grRateWT,...
    singleKOsim.mixo.hasEffect,...
    singleKOsim.mixo.delRxns,...
    singleKOsim.mixo.fluxSolution] = singleGeneDeletion(model_mixo, 'FBA', geneKOlist);

disp(singleKOsim.mixo.grRateWT);

for i = 1:length(singleKOsim.mixo.grRatio)
    if singleKOsim.mixo.grRatio(i) > 1.000001
        disp(singleKOsim.mixo.grRatio(i));
    end
end