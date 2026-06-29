function [model, appliedMedium, warnings] = applySingleModelMedium(model, mediumFile, options)
%APPLYSINGLEMODELMEDIUM Apply medium to a single-strain model.

if nargin < 3 || isempty(options)
    options = struct();
end
options = defaultMediumOptions(options);
warnings = {};
exchange = detectExchangeReactions(model);

if options.close_unlisted_uptakes
    for i = 1:numel(exchange.exchangeRxns)
        rxnIdx = find(strcmp(model.rxns, exchange.exchangeRxns{i}), 1);
        if ~isempty(rxnIdx) && model.lb(rxnIdx) < 0
            model.lb(rxnIdx) = 0;
        end
    end
end

medium = readMediumFileLocal(mediumFile);
applied = false(height(medium), 1);
for i = 1:height(medium)
    rxn = char(string(medium.exchange_rxn(i)));
    if isempty(rxn) || strcmpi(rxn, 'missing')
        rxn = findMetaboliteExchange(model, {char(string(medium.metabolite(i)))});
    end
    rxnIdx = find(strcmp(model.rxns, rxn), 1);
    if isempty(rxnIdx) && startsWith(rxn, 'EX_')
        rxnIdx = find(strcmp(model.rxns, ['R_' rxn]), 1);
        if ~isempty(rxnIdx)
            rxn = ['R_' rxn];
        end
    end
    if isempty(rxnIdx) && startsWith(rxn, 'R_EX_')
        strippedRxn = char(extractAfter(rxn, 2));
        rxnIdx = find(strcmp(model.rxns, strippedRxn), 1);
        if ~isempty(rxnIdx)
            rxn = strippedRxn;
        end
    end
    if isempty(rxnIdx)
        warnings{end+1} = sprintf('Single-model medium exchange reaction not found: %s', rxn); %#ok<AGROW>
        continue
    end
    lb = medium.lower_bound(i);
    ub = medium.upper_bound(i);
    if strcmpi(options.condition, 'anaerobic') && any(strcmpi(stripRPrefix(rxn), {'EX_o2_e','EX_o2(e)','EX_o2'}))
        lb = 0;
    elseif strcmpi(options.condition, 'microaerobic') && any(strcmpi(stripRPrefix(rxn), {'EX_o2_e','EX_o2(e)','EX_o2'}))
        lb = max(lb, -1);
    end
    model.lb(rxnIdx) = lb;
    model.ub(rxnIdx) = ub;
    applied(i) = true;
end

appliedMedium = medium(applied, :);
if ~isempty(appliedMedium)
    appliedMedium.applied = true(height(appliedMedium), 1);
end
end

function options = defaultMediumOptions(options)
if ~isfield(options, 'close_unlisted_uptakes')
    options.close_unlisted_uptakes = true;
end
if ~isfield(options, 'condition') || isempty(options.condition)
    options.condition = 'anaerobic';
end
end

function rxn = stripRPrefix(rxn)
rxn = regexprep(char(rxn), '^R_', '');
end

function medium = readMediumFileLocal(pathName)
[~, ~, ext] = fileparts(pathName);
if any(strcmpi(ext, {'.tsv', '.csv', '.txt'}))
    medium = readTsvOrCsv(pathName);
elseif strcmpi(ext, '.json')
    raw = jsondecode(fileread(pathName));
    medium = struct2table(raw);
elseif any(strcmpi(ext, {'.yml', '.yaml'}))
    error('applySingleModelMedium:YamlMediumUnsupported', 'Use TSV/CSV/JSON for medium files.');
else
    error('applySingleModelMedium:UnsupportedMedium', 'Unsupported medium file: %s', pathName);
end
required = {'metabolite','exchange_rxn','lower_bound','upper_bound'};
if ~all(ismember(required, medium.Properties.VariableNames))
    error('applySingleModelMedium:InvalidMedium', 'Medium table must contain metabolite, exchange_rxn, lower_bound, upper_bound.');
end
end
