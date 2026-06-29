function fluxes = extractDenitrificationFluxes(model, solution, aliasTable)
%EXTRACTDENITRIFICATIONFLUXES Convert exchange fluxes into positive metrics.
%
%   uptake = max(0, -exchangeFlux), secretion = max(0, exchangeFlux), and net
%   flux is the raw COBRA exchange flux.

if nargin < 3
    aliasTable = table();
end
targets = {'nitrate','nitrite','nitric_oxide','nitrous_oxide','dinitrogen','ammonium','oxygen','carbon_dioxide'};
canonical = {'no3','no2','no','n2o','n2','nh4','o2','co2'};
for i = 1:numel(targets)
    aliases = aliasesForTarget(aliasTable, canonical{i}, targets{i});
    rxn = findMetaboliteExchange(model, aliases);
    v = NaN;
    if ~isempty(rxn) && isfield(solution, 'x') && numel(solution.x) == numel(model.rxns)
        idx = find(strcmp(model.rxns, rxn), 1);
        v = solution.x(idx);
    end
    fluxes.(targets{i}).exchange_rxn = rxn;
    fluxes.(targets{i}).net_flux = v;
    fluxes.(targets{i}).uptake = max(0, -v);
    fluxes.(targets{i}).secretion = max(0, v);
end
end

function aliases = aliasesForTarget(aliasTable, canonical, category)
aliases = {canonical, [canonical '_e'], ['EX_' canonical '_e'], category};
if ~isempty(aliasTable) && all(ismember({'canonical_id','alias'}, aliasTable.Properties.VariableNames))
    rows = strcmp(string(aliasTable.canonical_id), canonical) | strcmp(string(aliasTable.category), category);
    aliases = unique([aliases, cellstr(string(aliasTable.alias(rows)))'], 'stable');
end
end
