function deleteCompiledModel(modelName)
models_dir = fullfile( fileparts(mfilename('fullpath')), 'stanModels');
modelFile = fullfile(models_dir, [modelName '.stan']);

delete(strrep(modelFile,'.stan','.exe'));
delete(strrep(modelFile,'.stan','.mat'));
delete(strrep(modelFile,'.stan','.hpp'));
end