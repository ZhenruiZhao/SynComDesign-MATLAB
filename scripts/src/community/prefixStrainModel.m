function result = prefixStrainModel(model, strainName, sharedCompartment)
%PREFIXSTRAINMODEL Prefix strain IDs and connect exchange metabolites to a shared pool.

safeName = matlab.lang.makeValidName(strainName);
[exchangeRxns, exchangeMets] = boundaryExchangeReactions(model);

newMets = cell(size(model.mets));
metRows = cell(numel(model.mets), 4);
for i = 1:numel(model.mets)
    sourceMet = char(model.mets{i});
    exIdx = find(strcmp(exchangeMets, sourceMet), 1);
    if ~isempty(exIdx)
        newMets{i} = [safeName, '__', sourceMet];
        role = 'strain_exchange_metabolite';
    else
        newMets{i} = [safeName, '__', sourceMet];
        role = 'strain_internal';
    end
    metRows(i, :) = {newMets{i}, strainName, sourceMet, role};
end

newRxns = cell(size(model.rxns));
rxnRows = cell(numel(model.rxns), 4);
for i = 1:numel(model.rxns)
    sourceRxn = char(model.rxns{i});
    newRxns{i} = [safeName, '__', sourceRxn];
    if any(strcmp(exchangeRxns, sourceRxn))
        role = 'strain_shared_interface';
    else
        role = 'strain_internal';
    end
    rxnRows(i, :) = {newRxns{i}, strainName, sourceRxn, role};
end

prefixed = model;
prefixed.mets = newMets(:);
prefixed.rxns = newRxns(:);
prefixed = addSharedInterfaceRows(prefixed, model, exchangeRxns, exchangeMets, sharedCompartment);

biomassRxn = '';
if isfield(model, 'syncomdesign') && isfield(model.syncomdesign, 'biomassRxn') && ~isempty(model.syncomdesign.biomassRxn)
    biomassRxn = [safeName, '__', model.syncomdesign.biomassRxn];
end

result = struct('model', prefixed, 'reactionMap', ...
    cell2table(rxnRows, 'VariableNames', {'community_rxn','strain','source_rxn','role'}), ...
    'metaboliteMap', cell2table(metRows, 'VariableNames', {'community_met','strain','source_met','role'}), ...
    'biomassRxn', biomassRxn);
end

function prefixed = addSharedInterfaceRows(prefixed, sourceModel, exchangeRxns, exchangeMets, sharedCompartment)
for i = 1:numel(exchangeRxns)
    sourceRxn = exchangeRxns{i};
    sourceRxnIdx = find(strcmp(sourceModel.rxns, sourceRxn), 1);
    if isempty(sourceRxnIdx)
        continue
    end
    sourceMet = exchangeMets{i};
    sourceMetIdx = find(strcmp(sourceModel.mets, sourceMet), 1);
    if isempty(sourceMetIdx)
        sharedCoeff = 1;
    else
        coeff = full(sourceModel.S(sourceMetIdx, sourceRxnIdx));
        if coeff == 0
            sharedCoeff = 1;
        else
            sharedCoeff = -coeff;
        end
    end
    sharedMet = sprintf('%s[%s]', canonicalMetId(sourceMet), sharedCompartment);
    sharedIdx = find(strcmp(prefixed.mets, sharedMet), 1);
    if isempty(sharedIdx)
        prefixed.S(end+1, :) = sparse(1, size(prefixed.S, 2));
        prefixed.mets{end+1, 1} = sharedMet;
        sharedIdx = numel(prefixed.mets);
    end
    prefixed.S(sharedIdx, sourceRxnIdx) = prefixed.S(sharedIdx, sourceRxnIdx) + sharedCoeff;
end
end

function [exchangeRxns, exchangeMets] = boundaryExchangeReactions(model)
rxns = cellstr(model.rxns);
isBoundary = (startsWith(rxns, 'EX_') | startsWith(rxns, 'R_EX_')) & ~isTransportReactionId(rxns);
idx = find(isBoundary);
exchangeRxns = rxns(idx);
exchangeMets = cell(numel(idx), 1);
for i = 1:numel(idx)
    metIdx = find(model.S(:, idx(i)) ~= 0, 1);
    if isempty(metIdx)
        exchangeMets{i} = exchangeMetFromRxnId(rxns{idx(i)});
    else
        exchangeMets{i} = char(model.mets{metIdx});
    end
end
valid = ~cellfun(@isempty, exchangeMets);
exchangeRxns = exchangeRxns(valid);
exchangeMets = exchangeMets(valid);
end

function met = exchangeMetFromRxnId(rxn)
met = regexprep(char(rxn), '^R_EX_', '');
met = regexprep(met, '^EX_', '');
met = regexprep(met, '_(e|u)$', '_e');
end

function tf = isTransportReactionId(rxns)
tf = false(size(rxns));
for i = 1:numel(rxns)
    tf(i) = ~isempty(regexp(char(rxns{i}), '(tex|tpp|t2pp|t3pp|tipp|tppi)$', 'once'));
end
end

function canonical = canonicalMetId(met)
canonical = regexprep(char(met), '\[[^\]]+\]$', '');
canonical = regexprep(canonical, '_[A-Za-z][A-Za-z0-9]*$', '');
canonical = regexprep(canonical, '^M_', '');
end
