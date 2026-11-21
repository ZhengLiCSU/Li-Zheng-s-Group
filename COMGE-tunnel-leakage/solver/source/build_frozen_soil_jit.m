function [kernel,info]=build_frozen_soil_jit(outputDir)
%BUILD_FROZEN_SOIL_JIT Build an opt-in JIT MEX without an external compiler.
% [kernel,info]=build_frozen_soil_jit() copies only the independent entry
% into an ASCII temporary folder and asks MATLAB Coder for a JIT MEX.
% The core and original public helpers are not changed or shadowed.
% Inspect the code-generation report: Build Type must be JIT MEX Function.
% MATLAB Coder and execution by the same MATLAB release are required.

assert(exist('codegen','file')~=0&&license('test','MATLAB_Coder'), ...
    'FrozenSoilJIT:CoderUnavailable','MATLAB Coder is required for JIT generation.');
if nargin<1||isempty(outputDir)
    outputDir=fullfile(tempdir,['frozen_soil_jit_' version('-release')]);
end
outputDir=char(outputDir);
if ~exist(outputDir,'dir'), mkdir(outputDir); end
% Resolve before changing the working directory.
originalDirectory=pwd;
cd(outputDir); outputDir=pwd; cd(originalDirectory);
assert(all(double(outputDir)<128),'FrozenSoilJIT:BuildPath', ...
    'Use an ASCII build folder; the default temp folder is suitable on this host.');
entry='frozen_soil_properties_jit'; binary='frozen_soil_properties_jit_mex';
source=fullfile(fileparts(mfilename('fullpath')),[entry '.m']);
copy=fullfile(outputDir,[entry '.m']);
assert(~strcmp(source,copy),'FrozenSoilJIT:BuildPath', ...
    'Choose a separate build folder, not the source folder.');
copyfile(source,copy);
restoreDirectory=onCleanup(@() cd(originalDirectory));
cd(outputDir);
vectorType=coder.typeof(0,[Inf,1],[true,false]);
% Keep the source/MEX folder separate from the removable code-generation
% folder. Request -jit explicitly with the documented CLI, rather than relying
% on a configuration-object setting whose report may precede final JIT build.
buildDirectory=fullfile(outputDir,'codegen_jit');
codegen(entry,'-config:mex','-jit', ...
    '-args',{vectorType,vectorType,zeros(22,1)}, ...
    '-o',binary,'-d',buildDirectory,'-report','-v');
clear restoreDirectory
addpath(outputDir,'-begin');
rehash;
expected=fullfile(outputDir,[binary '.' mexext]);
% Standard codegen saves the final MEX in the working folder, separately from
% -d. Resolve actual output explicitly if a release writes it inside -d.
if exist(expected,'file')==0
    generated=dir(fullfile(buildDirectory,'**',[binary '.' mexext]));
    assert(numel(generated)==1,'FrozenSoilJIT:BinaryResolution', ...
        ['JIT generation did not produce one callable MEX. Source/MEX folder: ', ...
         '%s; generated folder: %s. Inspect this build report and console.'], ...
        outputDir,buildDirectory);
    actualDirectory=generated(1).folder;
    addpath(actualDirectory,'-begin'); rehash;
    expected=fullfile(actualDirectory,generated(1).name);
end
assert(strcmp(which(binary),expected),'FrozenSoilJIT:BinaryResolution', ...
    'The compiled entry resolves to a different build.');
kernel=str2func(binary);
info=struct('requestedBackend','MATLAB Coder JIT MEX', ...
    'matlabRelease',version('-release'),'matlabVersion',version, ...
    'outputDir',outputDir,'buildDirectory',buildDirectory, ...
    'binaryFile',expected,'sourceFile',source, ...
    'sourceCopy',copy,'sourceText',fileread(source), ...
    'builtAtLocal',datestr(now,30),'entryName',entry,'binaryName',binary, ...
    'parameterCount',22,'propertyCount',12, ...
    'reportCheck','Verify report Summary Build Type is JIT MEX Function.');
save(fullfile(outputDir,'frozen_soil_jit_build_info.mat'),'info');
fprintf('Isolated JIT entry: %s\n',expected);
end
