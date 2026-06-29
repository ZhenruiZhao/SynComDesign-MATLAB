function [model, appliedMedium, warnings] = applyMedium(model, mediumFile, options)
%APPLYMEDIUM Legacy wrapper for single-model medium application.
%
% New code should call applySingleModelMedium for single-strain validation and
% applyCommunityExternalMedium for community FBA.

[model, appliedMedium, warnings] = applySingleModelMedium(model, mediumFile, options);
end
