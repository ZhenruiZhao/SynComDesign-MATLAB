function row = resultRowFromSolution(combinationId, strainNames, model, solution, objectiveMode, runtimeSeconds, aliasTable)
%RESULTROWFROMSOLUTION Build one stable community summary row.
if nargin < 7
    aliasTable = table();
end
totalBiomass = NaN;
strainBiomass = zeros(1, numel(strainNames));
if isfield(solution, 'f')
    totalBiomass = solution.f;
end
if isfield(model, 'syncomdesign') && isfield(model.syncomdesign, 'biomassMap') && isfield(solution, 'x')
    for i = 1:height(model.syncomdesign.biomassMap)
        idx = find(strcmp(model.rxns, char(model.syncomdesign.biomassMap.biomass_rxn(i))), 1);
        if ~isempty(idx) && idx <= numel(solution.x)
            strainBiomass(i) = solution.x(idx);
        end
    end
end
fluxes = extractDenitrificationFluxes(model, solution, aliasTable);
status = "unknown";
if isfield(solution, 'stat')
    status = string(solution.stat);
end
active = strainNames(strainBiomass > 1e-9);
row = table(string(combinationId), numel(strainNames), string(strjoin(strainNames, ';')), ...
    isFeasible(solution), status, string(objectiveMode), totalBiomass, ...
    string(strjoin(active, ';')), fluxes.nitrate.uptake, fluxes.nitrite.uptake, ...
    fluxes.nitrite.secretion, fluxes.nitric_oxide.uptake, fluxes.nitric_oxide.secretion, ...
    fluxes.nitrous_oxide.uptake, fluxes.nitrous_oxide.secretion, fluxes.nitrous_oxide.net_flux, ...
    fluxes.dinitrogen.secretion, safeDivide(fluxes.nitrate.uptake, totalBiomass), ...
    safeDivide(fluxes.nitrous_oxide.uptake, totalBiomass), safeDivide(fluxes.dinitrogen.secretion, totalBiomass), ...
    all(strainBiomass > 1e-9), runtimeSeconds, string(''), ...
    'VariableNames', {'combination_id','community_size','strain_names','feasible','solver_status', ...
    'objective_mode','total_biomass','active_strains','nitrate_uptake','nitrite_uptake', ...
    'nitrite_secretion','no_uptake','no_secretion','n2o_uptake','n2o_secretion','n2o_net_flux', ...
    'n2_secretion','nitrate_uptake_per_biomass','n2o_uptake_per_biomass','n2_production_per_biomass', ...
    'minimum_growth_satisfied','runtime_seconds','warning_message'});
for i = 1:numel(strainNames)
    row.(['strain_biomass_' matlab.lang.makeValidName(strainNames{i})]) = strainBiomass(i);
end
end

function ok = isFeasible(solution)
ok = isfield(solution, 'stat') && solution.stat == 1;
end
