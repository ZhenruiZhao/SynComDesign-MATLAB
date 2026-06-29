function summary = runSynComDesign(configPath)
%RUN_SYNCOMDESIGN_CARVEME Run SynComDesign on CarveMe/BiGG COBRA models.
%
%   summary = runSynComDesign('config/syncomdesign_config.yml') loads models,
%   applies a BiGG-style medium, enumerates strain combinations, builds shared
%   environment community models, solves configured objectives, and writes
%   result tables.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(rootDir));
if nargin < 1 || isempty(configPath)
    configPath = fullfile(rootDir, 'config', 'syncomdesign_config.yml');
end
config = loadSynComDesignConfig(configPath);
outputDir = fullfile(rootDir, char(config.project.output_dir));
ensureDirectory(outputDir);
logFile = fullfile(outputDir, 'run.log');
logMessage(logFile, 'Starting SynComDesign CarveMe/BiGG run.');
initializeCobraForSynComDesign(config, logFile);

modelDir = resolvePath(rootDir, config.models.directory);
files = dir(fullfile(modelDir, config.models.file_pattern));
models = struct('name', {}, 'path', {}, 'model', {}, 'biomassRxn', {}, 'exchange', {}, 'validation', {});
validationRows = {};
for i = 1:numel(files)
    pathName = fullfile(files(i).folder, files(i).name);
    [~, strainName] = fileparts(pathName);
    try
        info = loadAndStandardizeModel(pathName, config, strainName);
        models(end+1) = info; %#ok<AGROW>
        warnText = strjoin(info.validation.warnings, '; ');
        validationRows(end+1, :) = {strainName, pathName, info.validation.isValid, info.biomassRxn, warnText}; %#ok<AGROW>
    catch err
        validationRows(end+1, :) = {strainName, pathName, false, '', err.message}; %#ok<AGROW>
        logMessage(logFile, ['Model load failed: ' err.message]);
    end
end
validationTable = cell2table(validationRows, 'VariableNames', ...
    {'strain','model_path','valid','biomass_rxn','warning_message'});

species = 1:numel(models);
if isempty(config.combinations.max_size)
    maxSize = numel(species);
else
    maxSize = config.combinations.max_size;
end
requiredSpecies = normalizeSpeciesIndexList(config.combinations.required_species);
excludedSpecies = normalizeSpeciesIndexList(config.combinations.excluded_species);
if strcmpi(string(config.objective.type), "target_strain_biomass")
    targetIdx = find(strcmpi(string({models.name}), string(config.objective.target_strain)), 1);
    if isempty(targetIdx)
        error('runSynComDesign:UnknownTargetStrain', ...
            'target_strain "%s" was not found among loaded models: %s.', ...
            char(config.objective.target_strain), strjoin({models.name}, ', '));
    end
    requiredSpecies = unique([requiredSpecies(:); targetIdx], 'stable')';
end
[combinations, comboInfo] = GetAllCombination(species, ...
    'minCommunitySize', config.combinations.min_size, ...
    'maxCommunitySize', maxSize, ...
    'requiredSpecies', requiredSpecies, ...
    'excludedSpecies', excludedSpecies, ...
    'maxCombinations', config.combinations.max_combinations);
logMessage(logFile, sprintf('Total combinations: %d', comboInfo.totalCombinations));

aliasTable = mapMetaboliteAliases(resolvePath(rootDir, config.models.metabolite_aliases_file));
communityRows = table();
singleRows = table();
failedRows = table();
reactionMap = table();
metaboliteMap = table();

for i = 1:numel(combinations)
    combo = combinations{i};
    comboModels = models(combo);
    strainNames = {comboModels.name};
    comboId = strjoin(strainNames, '+');
    timer = tic;
    try
        community = buildCommunityModel(comboModels, config);
        mediumOptions = config.medium;
        mediumOptions.shared_environment_compartment = config.community.shared_environment_compartment;
        if strcmpi(string(config.medium.community_medium_mode), "legacy_all_exchange")
            [community, ~] = applyMedium(community, resolvePath(rootDir, config.medium.file), config.medium);
        else
            [community, ~] = applyCommunityExternalMedium(community, resolvePath(rootDir, config.medium.file), mediumOptions);
        end
        community = setCommunityObjective(community, config.objective);
        if any(strcmpi(string(config.objective.type), ["fixed_composition","equal_composition"]))
            ratios = [];
            if isfield(config.objective, 'composition_ratio')
                ratios = config.objective.composition_ratio;
            end
            if ~isempty(ratios) && numel(ratios) == numel(models)
                ratios = ratios(combo);
            end
            community = addFixedCompositionConstraint(community, ratios);
        end
        if config.community.require_all_species_active
            community = addAllSpeciesActiveConstraint(community, config.community.minimum_biomass_flux);
        end
        solution = solveConfiguredObjective(community, config.objective, aliasTable);
        row = resultRowFromSolution(comboId, strainNames, community, solution, config.objective.type, toc(timer), aliasTable);
        communityRows = appendTable(communityRows, row);
        reactionMap = appendTable(reactionMap, community.syncomdesign.reactionMap);
        metaboliteMap = appendTable(metaboliteMap, community.syncomdesign.metaboliteMap);
        if numel(combo) == 1
            singleRows = appendTable(singleRows, normalizeSingleRow(row, strainNames{1}));
        end
        if config.analysis.save_community_models
            save(fullfile(outputDir, ['community_' matlab.lang.makeValidName(comboId) '.mat']), 'community', 'solution');
        end
    catch err
        failed = table(string(comboId), string(err.message), toc(timer), ...
            'VariableNames', {'combination_id','error_message','runtime_seconds'});
        failedRows = appendTable(failedRows, failed);
        logMessage(logFile, ['Combination failed: ' comboId ' :: ' err.message]);
        if ~config.analysis.continue_on_error
            rethrow(err);
        end
    end
end

outputs = struct();
outputs.community_summary = addComparisonColumns(communityRows, singleRows);
outputs.single_strain_results = singleRows;
outputs.flux_ranges = table();
outputs.model_validation = validationTable;
outputs.reaction_mapping = unique(reactionMap, 'rows');
outputs.metabolite_mapping = unique(metaboliteMap, 'rows');
outputs.failed_combinations = failedRows;
writeSynComDesignOutputs(outputDir, outputs);

summary = makeRunSummary(models, outputs, comboInfo.totalCombinations, outputDir);
printRunSummary(summary);
logMessage(logFile, 'Finished SynComDesign CarveMe/BiGG run.');
end

function pathName = resolvePath(rootDir, pathName)
pathName = char(pathName);
if isempty(pathName)
    pathName = rootDir;
elseif ~isfolder(pathName) && ~isfile(pathName) && ~startsWith(pathName, rootDir) && isempty(regexp(pathName, '^[A-Za-z]:', 'once'))
    pathName = fullfile(rootDir, pathName);
end
end

function sol = solveModelOrFallback(model)
if exist('optimizeCbModel', 'file') == 2
    sol = optimizeCbModel(model);
else
    sol = struct('stat', 0, 'f', NaN, 'x', zeros(numel(model.rxns), 1));
end
end

function values = normalizeSpeciesIndexList(values)
if isempty(values)
    values = [];
elseif iscell(values)
    if isempty(values)
        values = [];
    elseif all(cellfun(@isnumeric, values))
        values = cell2mat(values);
    else
        values = str2double(string(values));
    end
elseif isstring(values) || ischar(values)
    values = str2double(string(values));
end
values = values(:)';
values = values(~isnan(values));
end

function sol = solveConfiguredObjective(model, objectiveConfig, aliasTable)
if strcmpi(string(objectiveConfig.type), "growth_then_n2o_consumption")
    sol = solveGrowthThenN2OUptake(model, objectiveConfig, aliasTable);
else
    sol = solveModelOrFallback(model);
end
end

function sol = solveGrowthThenN2OUptake(model, objectiveConfig, aliasTable)
growthModel = setCommunityObjective(model, struct('type', 'total_biomass'));
growthSol = solveModelOrFallback(growthModel);
if ~isFeasibleSolution(growthSol)
    sol = growthSol;
    return
end
maxGrowth = totalBiomassFlux(growthModel, growthSol);
if isnan(maxGrowth)
    maxGrowth = growthSol.f;
end
growthFraction = 0.9;
if isfield(objectiveConfig, 'growth_fraction') && ~isempty(objectiveConfig.growth_fraction)
    growthFraction = objectiveConfig.growth_fraction;
end
targetRxn = findN2OExchangeReaction(growthModel, aliasTable);
if isempty(targetRxn)
    sol = growthSol;
    sol.f = maxGrowth;
    return
end
functionalModel = addMinimumTotalBiomassConstraint(growthModel, growthFraction * maxGrowth);
functionalModel.c(:) = 0;
targetIdx = find(strcmp(functionalModel.rxns, targetRxn), 1);
functionalModel.c(targetIdx) = -1;
functionalSol = solveModelOrFallback(functionalModel);
growthN2OUptake = exchangeUptake(growthModel, growthSol, targetRxn);
functionalN2OUptake = exchangeUptake(functionalModel, functionalSol, targetRxn);
if isFeasibleSolution(functionalSol) && functionalN2OUptake > growthN2OUptake + 1e-9
    sol = functionalSol;
    sol.f = totalBiomassFlux(functionalModel, functionalSol);
    if isnan(sol.f)
        sol.f = maxGrowth;
    end
else
    sol = growthSol;
    sol.f = maxGrowth;
end
end

function uptake = exchangeUptake(model, sol, rxn)
uptake = 0;
if ~isFeasibleSolution(sol) || isempty(rxn)
    return
end
idx = find(strcmp(model.rxns, rxn), 1);
if ~isempty(idx) && idx <= numel(sol.x)
    uptake = max(0, -sol.x(idx));
end
end

function ok = isFeasibleSolution(sol)
ok = isfield(sol, 'stat') && sol.stat == 1 && isfield(sol, 'x') && ~isempty(sol.x);
end

function value = totalBiomassFlux(model, sol)
value = NaN;
if ~isfield(model, 'syncomdesign') || ~isfield(model.syncomdesign, 'biomassMap') || ~isfield(sol, 'x')
    return
end
value = 0;
biomassRxns = cellstr(model.syncomdesign.biomassMap.biomass_rxn);
for i = 1:numel(biomassRxns)
    idx = find(strcmp(model.rxns, biomassRxns{i}), 1);
    if ~isempty(idx) && idx <= numel(sol.x)
        value = value + sol.x(idx);
    end
end
end

function model = addMinimumTotalBiomassConstraint(model, minimumGrowth)
if isnan(minimumGrowth) || minimumGrowth <= 0
    return
end
biomassRxns = cellstr(model.syncomdesign.biomassMap.biomass_rxn);
biomassIdx = [];
for i = 1:numel(biomassRxns)
    idx = find(strcmp(model.rxns, biomassRxns{i}), 1);
    if ~isempty(idx)
        biomassIdx(end+1) = idx; %#ok<AGROW>
    end
end
if isempty(biomassIdx)
    return
end
constraintMet = 'syncomdesign_min_total_biomass_constraint';
constraintRxn = 'SYNCOMDESIGN_min_total_biomass_sink';
model.S(end+1, :) = sparse(1, size(model.S, 2));
model.S(end, biomassIdx) = 1;
newColumn = sparse(size(model.S, 1), 1);
newColumn(end) = -1;
model.S(:, end+1) = newColumn;
model.mets{end+1, 1} = constraintMet;
model.rxns{end+1, 1} = constraintRxn;
model.lb(end+1, 1) = minimumGrowth;
model.ub(end+1, 1) = max(1000, minimumGrowth * 10);
model.c(end+1, 1) = 0;
model = syncCobraConstraintFields(model);
end

function rxn = findN2OExchangeReaction(model, aliasTable)
aliases = {'n2o','n2o_e','EX_n2o_e','nitrous_oxide'};
if ~isempty(aliasTable) && all(ismember({'canonical_id','alias'}, aliasTable.Properties.VariableNames))
    rows = strcmp(string(aliasTable.canonical_id), 'n2o') | strcmp(string(aliasTable.category), 'nitrous_oxide');
    aliases = unique([aliases, cellstr(string(aliasTable.alias(rows)))'], 'stable');
end
rxn = findMetaboliteExchange(model, aliases);
end

function out = appendTable(out, row)
if isempty(out)
    out = row;
else
    [out, row] = alignTablesForAppend(out, row);
    out = [out; row]; %#ok<AGROW>
end
end

function [a, b] = alignTablesForAppend(a, b)
varsA = a.Properties.VariableNames;
varsB = b.Properties.VariableNames;
missingInA = setdiff(varsB, varsA, 'stable');
missingInB = setdiff(varsA, varsB, 'stable');
for i = 1:numel(missingInA)
    name = missingInA{i};
    a.(name) = missingColumnLike(b.(name), height(a));
end
for i = 1:numel(missingInB)
    name = missingInB{i};
    b.(name) = missingColumnLike(a.(name), height(b));
end
b = b(:, a.Properties.VariableNames);
end

function col = missingColumnLike(example, n)
if islogical(example)
    col = false(n, 1);
elseif isnumeric(example)
    col = NaN(n, 1);
elseif isstring(example)
    col = strings(n, 1);
elseif iscell(example)
    col = cell(n, 1);
else
    col = cell(n, 1);
end
end

function row = normalizeSingleRow(row, strain)
row.strain = string(strain);
end

function rows = addComparisonColumns(rows, singles)
if isempty(rows)
    return
end
defaults = {'best_member_biomass','expected_additive_biomass','biomass_gain_vs_best', ...
    'biomass_interaction_effect','n2o_gain_vs_best','nitrate_gain_vs_best'};
for d = 1:numel(defaults)
    rows.(defaults{d}) = NaN(height(rows), 1);
end
if isempty(singles) || ~ismember('strain', singles.Properties.VariableNames)
    return
end
for i = 1:height(rows)
    result = table2struct(rows(i, :));
    metrics = compareCommunityToSingles(result, singles);
    for d = 1:numel(defaults)
        rows.(defaults{d})(i) = metrics.(defaults{d});
    end
end
end

function summary = makeRunSummary(models, outputs, totalCombinations, outputDir)
summary = struct();
summary.ModelsDetected = numel(models);
summary.ModelsValidated = sum(outputs.model_validation.valid);
summary.BiomassReactionsDetected = sum(strlength(string(outputs.model_validation.biomass_rxn)) > 0);
summary.TotalCombinations = totalCombinations;
summary.FeasibleCombinations = countTrue(outputs.community_summary, 'feasible');
summary.FailedCombinations = height(outputs.failed_combinations);
summary.BestBiomassCombination = bestBy(outputs.community_summary, 'total_biomass');
summary.BestNitrateConsumingCombination = bestBy(outputs.community_summary, 'nitrate_uptake');
summary.BestN2OConsumingCombination = bestBy(outputs.community_summary, 'n2o_uptake');
summary.BestN2ProducingCombination = bestBy(outputs.community_summary, 'n2_secretion');
summary.OutputDirectory = outputDir;
end

function n = countTrue(tbl, name)
if isempty(tbl) || ~ismember(name, tbl.Properties.VariableNames)
    n = 0;
else
    n = sum(tbl.(name));
end
end

function combo = bestBy(tbl, metric)
combo = "";
if isempty(tbl) || ~ismember(metric, tbl.Properties.VariableNames)
    return
end
[value, idx] = max(tbl.(metric));
if ~isempty(idx) && ~isnan(value)
    combo = tbl.combination_id(idx);
end
end

function printRunSummary(summary)
fprintf('Models detected: %d\n', summary.ModelsDetected);
fprintf('Models validated: %d\n', summary.ModelsValidated);
fprintf('Biomass reactions detected: %d\n', summary.BiomassReactionsDetected);
fprintf('Total combinations: %d\n', summary.TotalCombinations);
fprintf('Feasible combinations: %d\n', summary.FeasibleCombinations);
fprintf('Failed combinations: %d\n', summary.FailedCombinations);
fprintf('Best biomass combination: %s\n', summary.BestBiomassCombination);
fprintf('Best nitrate-consuming combination: %s\n', summary.BestNitrateConsumingCombination);
fprintf('Best N2O-consuming combination: %s\n', summary.BestN2OConsumingCombination);
fprintf('Best N2-producing combination: %s\n', summary.BestN2ProducingCombination);
fprintf('Output directory: %s\n', summary.OutputDirectory);
end

function logMessage(logFile, message)
fid = fopen(logFile, 'a');
if fid > 0
    fprintf(fid, '[%s] %s\n', datestr(now, 31), message);
    fclose(fid);
end
end

function initializeCobraForSynComDesign(config, logFile)
if exist('initCobraToolbox', 'file') ~= 2
    logMessage(logFile, 'COBRA Toolbox initCobraToolbox not found on MATLAB path.');
    return
end
try
    initCobraToolbox(false);
    logMessage(logFile, 'COBRA Toolbox initialized.');
catch err
    logMessage(logFile, ['COBRA Toolbox initialization warning: ' err.message]);
end
if isfield(config, 'solver') && isfield(config.solver, 'name') && ~isempty(config.solver.name) ...
        && exist('changeCobraSolver', 'file') == 2
    solverName = char(config.solver.name);
    if any(strcmpi(solverName, {'auto','LP'}))
        solverName = '';
    end
    try
        if isempty(solverName)
            ok = false;
        else
            ok = changeCobraSolver(solverName, 'LP');
            logMessage(logFile, sprintf('Requested COBRA LP solver: %s (status=%g).', solverName, ok));
        end
    catch err
        logMessage(logFile, ['Requested solver failed: ' solverName ' :: ' err.message]);
        ok = false;
    end
    if ~ok
        fallbackSolvers = {'gurobi','ibm_cplex','mosek','glpk'};
        for i = 1:numel(fallbackSolvers)
            try
                ok = changeCobraSolver(fallbackSolvers{i}, 'LP');
                if ok
                    logMessage(logFile, ['Using fallback COBRA LP solver: ' fallbackSolvers{i}]);
                    return
                end
            catch
            end
        end
    end
end
end
