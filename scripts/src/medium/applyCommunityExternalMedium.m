function [model, appliedMedium, warnings, mapping] = applyCommunityExternalMedium(model, mediumFile, options)
%APPLYCOMMUNITYEXTERNALMEDIUM Apply medium only to external shared exchange.
%
% Closure targets come only from buildCommunityModel's externalExchangeMap.

if nargin < 3 || isempty(options)
    options = struct();
end
options = defaultCommunityMediumOptions(options);
warnings = {};

[externalRxns, externalMets] = externalSharedExchangeList(model);
if options.close_unlisted_external_medium_uptakes
    for i = 1:numel(externalRxns)
        rxnIdx = find(strcmp(model.rxns, externalRxns{i}), 1);
        if ~isempty(rxnIdx) && model.lb(rxnIdx) < 0
            model.lb(rxnIdx) = 0;
        end
    end
end

medium = readMediumFileLocal(mediumFile);
mappingRows = cell(height(medium), 8);
applied = false(height(medium), 1);
for i = 1:height(medium)
    mediumRxn = char(string(medium.exchange_rxn(i)));
    mediumMet = char(string(medium.metabolite(i)));
    sharedMet = expectedSharedMet(mediumRxn, mediumMet, options.shared_environment_compartment);
    mapIdx = find(strcmp(externalMets, sharedMet), 1);
    found = ~isempty(mapIdx);
    sharedRxn = '';
    if found
        sharedRxn = externalRxns{mapIdx};
    end

    lb = medium.lower_bound(i);
    ub = medium.upper_bound(i);
    if strcmpi(options.condition, 'anaerobic') && isOxygenMedium(mediumRxn, mediumMet)
        lb = 0;
    elseif strcmpi(options.condition, 'microaerobic') && isOxygenMedium(mediumRxn, mediumMet)
        lb = max(lb, -1);
    end

    warning = '';
    if found
        rxnIdx = find(strcmp(model.rxns, sharedRxn), 1);
        model.lb(rxnIdx) = lb;
        model.ub(rxnIdx) = ub;
        applied(i) = true;
    else
        warning = sprintf('No external shared exchange found for %s/%s; expected shared metabolite %s.', mediumRxn, mediumMet, sharedMet);
        warnings{end+1} = warning; %#ok<AGROW>
    end
    mappingRows(i, :) = {mediumRxn, mediumMet, sharedMet, sharedRxn, found, lb, ub, warning};
end

mapping = cell2table(mappingRows, 'VariableNames', ...
    {'medium_exchange_rxn','medium_metabolite','shared_metabolite','shared_external_exchange','found','lower_bound','upper_bound','warning'});
appliedMedium = medium(applied, :);
if ~isempty(appliedMedium)
    appliedMedium.applied = true(height(appliedMedium), 1);
end

model.syncomdesign.communityMedium = struct( ...
    'mode', options.community_medium_mode, ...
    'allow_cross_feeding', options.allow_cross_feeding, ...
    'close_unlisted_external_medium_uptakes', options.close_unlisted_external_medium_uptakes, ...
    'close_strain_interface_uptakes', options.close_strain_interface_uptakes, ...
    'close_internal_transport', options.close_internal_transport, ...
    'mapping', mapping);
end

function [externalRxns, externalMets] = externalSharedExchangeList(model)
if ~isfield(model, 'syncomdesign') || ~isfield(model.syncomdesign, 'externalExchangeMap')
    error('applyCommunityExternalMedium:MissingExternalExchangeMap', ...
        'Community model must contain syncomdesign.externalExchangeMap from buildCommunityModel.');
end
externalRxns = cellstr(string(model.syncomdesign.externalExchangeMap.external_exchange_rxn));
externalMets = cellstr(string(model.syncomdesign.externalExchangeMap.shared_metabolite));
end

function options = defaultCommunityMediumOptions(options)
options = defaultField(options, 'community_medium_mode', 'external_shared_only');
options = defaultField(options, 'close_unlisted_external_medium_uptakes', true);
options = defaultField(options, 'allow_cross_feeding', true);
options = defaultField(options, 'close_strain_interface_uptakes', false);
options = defaultField(options, 'close_internal_transport', false);
options = defaultField(options, 'legacy_close_unlisted_uptakes', false);
options = defaultField(options, 'condition', 'anaerobic');
options = defaultField(options, 'shared_environment_compartment', 'u');
end

function s = defaultField(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end

function medium = readMediumFileLocal(pathName)
[~, ~, ext] = fileparts(pathName);
if any(strcmpi(ext, {'.tsv', '.csv', '.txt'}))
    medium = readTsvOrCsv(pathName);
elseif strcmpi(ext, '.json')
    raw = jsondecode(fileread(pathName));
    medium = struct2table(raw);
else
    error('applyCommunityExternalMedium:UnsupportedMedium', 'Unsupported medium file: %s', pathName);
end
required = {'metabolite','exchange_rxn','lower_bound','upper_bound'};
if ~all(ismember(required, medium.Properties.VariableNames))
    error('applyCommunityExternalMedium:InvalidMedium', 'Medium table must contain metabolite, exchange_rxn, lower_bound, upper_bound.');
end
end

function sharedMet = expectedSharedMet(mediumRxn, mediumMet, sharedCompartment)
base = char(mediumMet);
if isempty(base) || strcmpi(base, 'missing')
    base = regexprep(char(mediumRxn), '^R_EX_', '');
    base = regexprep(base, '^EX_', '');
    base = regexprep(base, '_(e|u)$', '');
end
base = regexprep(base, '^M_', '');
base = regexprep(base, '\[[^\]]+\]$', '');
base = regexprep(base, '_(e|u)$', '');
sharedMet = sprintf('%s[%s]', base, sharedCompartment);
end

function tf = isOxygenMedium(mediumRxn, mediumMet)
rxn = regexprep(char(mediumRxn), '^R_', '');
met = regexprep(char(mediumMet), '_(e|u)$', '');
tf = any(strcmpi(rxn, {'EX_o2_e','EX_o2_u','EX_o2','o2'})) || any(strcmpi(met, {'EX_o2_e','EX_o2_u','EX_o2','o2'}));
end
