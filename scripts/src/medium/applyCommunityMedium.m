function [model, appliedMedium, warnings, mapping] = applyCommunityMedium(model, mediumFile, options)
%APPLYCOMMUNITYMEDIUM Compatibility wrapper.
%
% Community FBA should use applyCommunityExternalMedium, whose closure targets
% come only from buildCommunityModel's external shared exchange list.

[model, appliedMedium, warnings, mapping] = applyCommunityExternalMedium(model, mediumFile, options);
end
