function sm = getCompiledModel(modelName)
models_dir = fullfile( fileparts(mfilename('fullpath')));
modelFile = fullfile(models_dir, [modelName '.mat']);

if exist( modelFile, 'file')
    load(modelFile, 'sm');
    fprintf('Loaded model %s\n', modelName);
    
else %If it doesn't exist, compile and save it
    fprintf('Compiled model %s not found, compiling now...\n', modelName);
    
    %Get stan model file
    stanFile = strrep(modelFile, '.mat', '.stan');
    
    %Compile model
    sm = StanModel('file',stanFile,'working_dir',tempdir);
    sm.compile();
    
    %Save copy of compiled model .mat object
    save(modelFile,'sm');
end
end