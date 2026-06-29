function classification = classifyCommunityReactions(model)
%CLASSIFYCOMMUNITYREACTIONS Classify community reactions for medium handling.

rxns = cellstr(model.rxns);
classes = repmat({'unknown'}, numel(rxns), 1);
mediumApplies = false(numel(rxns), 1);
roleMap = makeRoleMap(model);

for i = 1:numel(rxns)
    rxn = rxns{i};
    role = '';
    if isKey(roleMap, rxn)
        role = roleMap(rxn);
    end

    if strcmp(role, 'external_medium_exchange') || isSharedExternalExchange(model, i)
        classes{i} = 'external_medium_exchange';
        mediumApplies(i) = true;
    elseif isTransportReaction(rxn)
        classes{i} = 'internal_transport';
    elseif strcmp(role, 'strain_shared_interface') || contains(rxn, '__R_EX_') || contains(rxn, '__EX_')
        classes{i} = 'strain_shared_interface';
    elseif strcmp(role, 'strain_internal') || strcmp(role, 'exchange')
        classes{i} = 'metabolic_reaction';
    elseif nnz(model.S(:, i)) > 1
        classes{i} = 'metabolic_reaction';
    end
end

classification = table(rxns(:), classes(:), mediumApplies(:), ...
    'VariableNames', {'reaction_id','classification','medium_applies'});
end

function roleMap = makeRoleMap(model)
roleMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
if ~isfield(model, 'syncomdesign') || ~isfield(model.syncomdesign, 'reactionMap')
    return
end
mapping = model.syncomdesign.reactionMap;
for i = 1:height(mapping)
    roleMap(char(mapping.community_rxn{i})) = char(mapping.role{i});
end
end

function tf = isSharedExternalExchange(model, idx)
rxn = char(model.rxns{idx});
tf = startsWith(rxn, 'R_EX_') && nnz(model.S(:, idx)) == 1;
end

function tf = isTransportReaction(rxn)
rxn = char(rxn);
tf = ~isempty(regexp(rxn, '(tex|tpp|t2pp|t3pp|tipp|tppi)$', 'once'));
end
