function validation = validateCarveMeModel(model, biomass, exchange, config)
%VALIDATECARVEMEMODEL Run structural and biological prechecks.

warnings = {};
errors = {};
[metCount, rxnCount] = size(model.S);
if numel(model.mets) ~= metCount
    errors{end+1} = 'S row count does not match model.mets.'; %#ok<AGROW>
end
if numel(model.rxns) ~= rxnCount || numel(model.lb) ~= rxnCount || numel(model.ub) ~= rxnCount
    errors{end+1} = 'S column count does not match rxns/lb/ub.'; %#ok<AGROW>
end
if isfield(model, 'c') && numel(model.c) ~= rxnCount
    errors{end+1} = 'Objective vector c does not match reaction count.'; %#ok<AGROW>
end
if any(model.lb > model.ub)
    errors{end+1} = 'One or more reactions have lower bound greater than upper bound.'; %#ok<AGROW>
end
if isempty(biomass.biomassRxn)
    warnings{end+1} = 'Biomass reaction is missing or ambiguous; configure biomass_reactions.tsv.'; %#ok<AGROW>
end
if isempty(exchange.exchangeRxns)
    warnings{end+1} = 'No exchange reactions detected.'; %#ok<AGROW>
end
if ~isfield(model, 'metFormulas')
    warnings{end+1} = 'metFormulas field is missing.'; %#ok<AGROW>
end
if ~isfield(model, 'metCharges')
    warnings{end+1} = 'metCharges field is missing.'; %#ok<AGROW>
end

targets = {'no3','no2','no','n2o','n2'};
missingTargets = {};
for i = 1:numel(targets)
    if isempty(findMetaboliteExchange(model, targetAliases(targets{i})))
        missingTargets{end+1} = targets{i}; %#ok<AGROW>
    end
end
if ~isempty(missingTargets)
    warnings{end+1} = ['Missing denitrification exchange targets: ', strjoin(missingTargets, ', ')]; %#ok<AGROW>
end

growth = NaN;
solverStatus = 'not_run';
if ~isempty(biomass.biomassRxn) && exist('optimizeCbModel', 'file') == 2
    try
        testModel = applySingleModelMedium(model, config.medium.file, config.medium);
        testModel.c(:) = 0;
        rxnIdx = find(strcmp(testModel.rxns, biomass.biomassRxn), 1);
        testModel.c(rxnIdx) = 1;
        sol = optimizeCbModel(testModel);
        if isfield(sol, 'f')
            growth = sol.f;
        end
        if isfield(sol, 'stat')
            solverStatus = num2str(sol.stat);
        end
    catch err
        solverStatus = ['error: ', err.message];
    end
end

validation = struct('isValid', isempty(errors), 'errors', {errors}, ...
    'warnings', {warnings}, 'growthOnConfiguredMedium', growth, ...
    'solverStatus', solverStatus);
end

function aliases = targetAliases(name)
switch lower(name)
    case 'no3'
        aliases = {'no3','no3_e','EX_no3_e','nitrate'};
    case 'no2'
        aliases = {'no2','no2_e','EX_no2_e','nitrite'};
    case 'no'
        aliases = {'no','no_e','EX_no_e','nitric_oxide'};
    case 'n2o'
        aliases = {'n2o','n2o_e','EX_n2o_e','nitrous_oxide'};
    case 'n2'
        aliases = {'n2','n2_e','EX_n2_e','dinitrogen'};
end
end
