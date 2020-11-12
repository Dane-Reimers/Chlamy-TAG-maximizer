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
model_mixo = readCbModel('../data/models/iCre1355_auto.xml');
load ../data/transcript_gene_map

objectives = {'PLDAGAT1819Z1819Z1819Z4', 'PLDAGAT1819Z1819Z1819Z3', 'PLDAGAT1819Z1819Z1819Z2', 'PLDAGAT1819Z1819Z1819Z1', 'PLDAGAT1819Z1819Z18111Z3'};
model_mixo = changeObjective(model_mixo, "Biomass_Chlamy_auto", 1);

% process list of KO candidate genes via transcript_gene_map
tmp = unique(transcript_gene_map(2:end,2));
geneKOlist = cell(length(tmp),1);
for i=1:length(tmp)
    geneKOlist{i} = transcript_gene_map(...
        strcmp(tmp(i,1),transcript_gene_map(:,2)),1);
end
clear tmp


[singleKOsim.mixo.grRatio,...
    singleKOsim.mixo.grRateKO,...
    singleKOsim.mixo.grRateWT,...
    singleKOsim.mixo.hasEffect,...
    singleKOsim.mixo.delRxns,...
    singleKOsim.mixo.fluxSolution] = singleGeneDeletion(model_mixo, 'FBA', geneKOlist);

for i = 1:length(singleKOsim.mixo.grRatio)
    if singleKOsim.mixo.grRatio(i) > 1.000001
        disp(singleKOsim.mixo.grRatio(i));
    end
end