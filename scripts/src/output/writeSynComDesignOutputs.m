function writeSynComDesignOutputs(outputDir, outputs)
%WRITESYNCOMDESIGNOUTPUTS Write SynComDesign result tables.
ensureDirectory(outputDir);
writeIfPresent(outputs, 'community_summary', outputDir);
writeIfPresent(outputs, 'single_strain_results', outputDir);
writeIfPresent(outputs, 'flux_ranges', outputDir);
writeIfPresent(outputs, 'model_validation', outputDir);
writeIfPresent(outputs, 'reaction_mapping', outputDir);
writeIfPresent(outputs, 'metabolite_mapping', outputDir);
writeIfPresent(outputs, 'community_medium_requirements', outputDir);
writeIfPresent(outputs, 'failed_combinations', outputDir);
end

function writeIfPresent(outputs, name, outputDir)
if isfield(outputs, name) && istable(outputs.(name))
    writetable(outputs.(name), fullfile(outputDir, [name '.tsv']), 'FileType', 'text', 'Delimiter', '\t');
    writetable(outputs.(name), fullfile(outputDir, [name '.csv']));
else
    empty = table();
    writetable(empty, fullfile(outputDir, [name '.tsv']), 'FileType', 'text', 'Delimiter', '\t');
end
end
