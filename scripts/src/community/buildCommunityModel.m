function community = buildCommunityModel(modelInfos, config)
%BUILDCOMMUNITYMODEL Build a prefixed shared-environment community model.
%
%   community = buildCommunityModel(modelInfos, config) combines the selected
%   strain models without using any strain as a base model. Strain exchange
%   reactions connect strain extracellular metabolites to a shared environment,
%   and explicit shared external exchanges are added for medium control.

if isempty(modelInfos)
    error('buildCommunityModel:NoModels', 'At least one model is required.');
end
sharedCompartment = config.community.shared_environment_compartment;
parts = cell(numel(modelInfos), 1);
mappingRows = {};
metRows = {};
biomassRows = {};
for i = 1:numel(modelInfos)
    parts{i} = prefixStrainModel(modelInfos(i).model, modelInfos(i).name, sharedCompartment);
    for r = 1:height(parts{i}.reactionMap)
        mappingRows(end+1, :) = table2cell(parts{i}.reactionMap(r, :)); %#ok<AGROW>
    end
    for m = 1:height(parts{i}.metaboliteMap)
        metRows(end+1, :) = table2cell(parts{i}.metaboliteMap(m, :)); %#ok<AGROW>
    end
    biomassRows(end+1, :) = {modelInfos(i).name, parts{i}.biomassRxn}; %#ok<AGROW>
end

allMets = {};
allRxns = {};
allLb = [];
allUb = [];
allC = [];
S = sparse(0, 0);
for i = 1:numel(parts)
    part = parts{i}.model;
    [allMets, S] = mergeStoichiometry(S, allMets, part.S, part.mets);
    allRxns = [allRxns; part.rxns(:)]; %#ok<AGROW>
    allLb = [allLb; part.lb(:)]; %#ok<AGROW>
    allUb = [allUb; part.ub(:)]; %#ok<AGROW>
    if isfield(part, 'c')
        allC = [allC; part.c(:)]; %#ok<AGROW>
    else
        allC = [allC; zeros(numel(part.rxns), 1)]; %#ok<AGROW>
    end
end

community = struct();
community.S = S;
community.mets = allMets(:);
community.rxns = allRxns(:);
community.lb = allLb;
community.ub = allUb;
community.c = allC;
community.description = 'SynComDesign CarveMe/BiGG community model';
community.syncomdesign = struct();
community.syncomdesign.strainNames = {modelInfos.name};
community.syncomdesign.reactionMap = cell2table(mappingRows, 'VariableNames', ...
    {'community_rxn','strain','source_rxn','role'});
community.syncomdesign.metaboliteMap = cell2table(metRows, 'VariableNames', ...
    {'community_met','strain','source_met','role'});
community.syncomdesign.biomassMap = cell2table(biomassRows, 'VariableNames', ...
    {'strain','biomass_rxn'});
community = addSharedExternalExchanges(community, sharedCompartment);
community.syncomdesign.transportMap = community.syncomdesign.reactionMap(strcmp(community.syncomdesign.reactionMap.role, 'strain_shared_interface'), :);
community = syncCobraConstraintFields(community);
community = validateCommunityModel(community);
end

function community = addSharedExternalExchanges(community, sharedCompartment)
sharedMets = community.mets(endsWith(community.mets, ['[' sharedCompartment ']']));
sharedMets = unique(sharedMets, 'stable');
rows = cell(numel(sharedMets), 4);
for i = 1:numel(sharedMets)
    sharedMet = sharedMets{i};
    rxn = sharedExchangeRxnId(sharedMet, sharedCompartment);
    if ~any(strcmp(community.rxns, rxn))
        metIdx = find(strcmp(community.mets, sharedMet), 1);
        community.S(:, end+1) = sparse(size(community.S, 1), 1);
        community.S(metIdx, end) = -1;
        community.rxns{end+1, 1} = rxn;
        community.lb(end+1, 1) = -1000;
        community.ub(end+1, 1) = 1000;
        community.c(end+1, 1) = 0;
    end
    rows(i, :) = {rxn, 'external', rxn, 'external_medium_exchange'};
end
externalMap = cell2table(rows, 'VariableNames', {'community_rxn','strain','source_rxn','role'});
community.syncomdesign.reactionMap = [community.syncomdesign.reactionMap; externalMap];
community.syncomdesign.externalExchangeMap = table(sharedMets(:), externalMap.community_rxn(:), ...
    'VariableNames', {'shared_metabolite','external_exchange_rxn'});
community.syncomdesign.externalSharedExchangeRxns = externalMap.community_rxn(:);
end

function rxn = sharedExchangeRxnId(sharedMet, sharedCompartment)
base = regexprep(char(sharedMet), ['\[' sharedCompartment '\]$'], '');
base = regexprep(base, '[^A-Za-z0-9_]', '_');
rxn = sprintf('R_EX_%s_%s', base, sharedCompartment);
end

function [allMets, Smerged] = mergeStoichiometry(Smerged, allMets, Spart, metsPart)
oldRxnCount = size(Smerged, 2);
for j = 1:numel(metsPart)
    if ~ismember(metsPart{j}, allMets)
        allMets{end+1, 1} = metsPart{j}; %#ok<AGROW>
        Smerged(end+1, :) = sparse(1, size(Smerged, 2));
    end
end
newCols = sparse(numel(allMets), size(Spart, 2));
for j = 1:numel(metsPart)
    targetIdx = find(strcmp(allMets, metsPart{j}), 1);
    newCols(targetIdx, :) = Spart(j, :);
end
Smerged(:, oldRxnCount+1:oldRxnCount+size(Spart, 2)) = newCols;
end
