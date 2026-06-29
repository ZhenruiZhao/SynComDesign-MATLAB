function config = loadSynComDesignConfig(configPath)
%LOADSYNCOMDESIGNCONFIG Load SynComDesign JSON or simple YAML configuration.
%
%   config = loadSynComDesignConfig(configPath) returns a struct. JSON is parsed
%   with jsondecode. YAML support covers the simple nested key/value/list shape
%   used by config/syncomdesign_config.yml.

if nargin < 1 || isempty(configPath)
    configPath = fullfile('config', 'syncomdesign_config.json');
end
if ~isfile(configPath)
    error('loadSynComDesignConfig:MissingFile', 'Configuration file not found: %s', configPath);
end
[~, ~, ext] = fileparts(configPath);
text = fileread(configPath);
if ~isempty(text) && double(text(1)) == 65279
    text = text(2:end);
end
if strcmpi(ext, '.json')
    config = jsondecode(text);
elseif any(strcmpi(ext, {'.yml', '.yaml'}))
    config = parseSimpleYaml(text);
else
    error('loadSynComDesignConfig:UnsupportedFormat', 'Use .json, .yml, or .yaml configuration files.');
end
config = applyConfigDefaults(config);
end

function config = parseSimpleYaml(text)
lines = regexp(text, '\r?\n', 'split');
config = struct();
section = '';
for i = 1:numel(lines)
    line = strtrim(regexprep(lines{i}, '#.*$', ''));
    if isempty(line)
        continue
    end
    if endsWith(line, ':')
        section = matlab.lang.makeValidName(line(1:end-1));
        config.(section) = struct();
        continue
    end
    parts = regexp(line, '^([A-Za-z0-9_]+):\s*(.*)$', 'tokens', 'once');
    if isempty(parts)
        continue
    end
    key = matlab.lang.makeValidName(parts{1});
    value = parseYamlValue(strtrim(parts{2}));
    if isempty(section)
        config.(key) = value;
    else
        config.(section).(key) = value;
    end
end
end

function value = parseYamlValue(raw)
if isempty(raw) || strcmpi(raw, 'null')
    value = [];
elseif strcmpi(raw, 'true')
    value = true;
elseif strcmpi(raw, 'false')
    value = false;
elseif startsWith(raw, '[') && endsWith(raw, ']')
    body = strtrim(raw(2:end-1));
    if isempty(body)
        value = {};
        return
    end
    items = strtrim(strsplit(body, ','));
    numeric = str2double(items);
    if all(~isnan(numeric))
        value = numeric;
    else
        value = cellfun(@stripQuotes, items, 'UniformOutput', false);
    end
else
    number = str2double(raw);
    if ~isnan(number)
        value = number;
    else
        value = stripQuotes(raw);
    end
end
end

function out = stripQuotes(in)
out = char(in);
if numel(out) >= 2
    if (out(1) == '"' && out(end) == '"') || (out(1) == '''' && out(end) == '''')
        out = out(2:end-1);
    end
end
end

function config = applyConfigDefaults(config)
config = ensureSection(config, 'project');
config.project = defaultField(config.project, 'name', 'denitrification_community');
config.project = defaultField(config.project, 'output_dir', 'results');

config = ensureSection(config, 'models');
config.models = defaultField(config.models, 'directory', 'models');
config.models = defaultField(config.models, 'file_pattern', '*.xml');
config.models = defaultField(config.models, 'biomass_reactions_file', fullfile('config', 'biomass_reactions.tsv'));
config.models = defaultField(config.models, 'metabolite_aliases_file', fullfile('config', 'metabolite_aliases.tsv'));

config = ensureSection(config, 'combinations');
config.combinations = defaultField(config.combinations, 'min_size', 1);
config.combinations = defaultField(config.combinations, 'max_size', []);
config.combinations = defaultField(config.combinations, 'required_species', {});
config.combinations = defaultField(config.combinations, 'excluded_species', {});
config.combinations = defaultField(config.combinations, 'max_combinations', 100000);

config = ensureSection(config, 'medium');
config.medium = defaultField(config.medium, 'file', fullfile('config', 'medium.tsv'));
config.medium = defaultField(config.medium, 'condition', 'anaerobic');
config.medium = defaultField(config.medium, 'close_unlisted_uptakes', true);
config.medium = defaultField(config.medium, 'community_medium_mode', 'external_shared_only');
config.medium = defaultField(config.medium, 'close_unlisted_external_medium_uptakes', true);
config.medium = defaultField(config.medium, 'allow_cross_feeding', true);
config.medium = defaultField(config.medium, 'close_strain_interface_uptakes', false);
config.medium = defaultField(config.medium, 'close_internal_transport', false);
config.medium = defaultField(config.medium, 'legacy_close_unlisted_uptakes', false);

config = ensureSection(config, 'community');
config.community = defaultField(config.community, 'require_all_species_active', true);
config.community = defaultField(config.community, 'minimum_biomass_flux', 1e-6);
config.community = defaultField(config.community, 'shared_environment_compartment', 'u');

config = ensureSection(config, 'objective');
config.objective = defaultField(config.objective, 'type', 'growth_then_n2o_consumption');
config.objective = applyScenarioId(config.objective);
config.objective = defaultField(config.objective, 'growth_fraction', 0.9);
config.objective = defaultField(config.objective, 'target_strain', []);
config.objective = defaultField(config.objective, 'biomass_weights', 'equal');

config = ensureSection(config, 'analysis');
config.analysis = defaultField(config.analysis, 'run_fva', true);
config.analysis = defaultField(config.analysis, 'fva_fraction_of_optimum', 90);
config.analysis = defaultField(config.analysis, 'save_community_models', false);
config.analysis = defaultField(config.analysis, 'continue_on_error', true);
end

function objective = applyScenarioId(objective)
if ~isfield(objective, 'scenario_id') || isempty(objective.scenario_id)
    return
end
switch objective.scenario_id
    case 1
        objective.type = 'total_biomass';
    case 2
        objective.type = 'target_strain_biomass';
    case 3
        objective.type = 'equal_composition';
    case 4
        objective.type = 'fixed_composition';
    case 5
        objective.type = 'growth_then_n2o_consumption';
    otherwise
        error('loadSynComDesignConfig:UnknownScenarioId', 'objective.scenario_id must be 1, 2, 3, 4, or 5.');
end
end

function s = ensureSection(s, name)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = struct();
end
end

function s = defaultField(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
