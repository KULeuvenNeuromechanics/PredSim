function [pathBuildExternalFunction] = buildExternalFunction(filename, cpp_dir, compiler, pathMain)
%buildExternalFunction Builds external function based on opensim AD
%modified version of the opensims source code using the installed
%libraries.
%   filename: name of the cpp file
%   cpp_dir: folder with cpp file
%   compiler: e.g. Visual Studio 15 2017 Win64
%   pathMain: path to the main directory of the opensimAD folder (https://github.com/Lars-DHondt-KUL/opensimAD)

%% Part 1: create executable from .cpp file

pathBuildExpressionGraph = fullfile(pathMain, 'buildExpressionGraph');
pathBuild = fullfile(pathBuildExpressionGraph,filename);
if ~isfolder(pathBuild)
    mkdir(pathBuild)
end

OpenSimAD_DIR = fullfile(pathMain, 'OpenSimAD-install');
SDK_DIR = fullfile(OpenSimAD_DIR, 'sdk');
BIN_DIR = fullfile(OpenSimAD_DIR, 'bin');

% run cmake to create expression graph
cd(pathBuild)
cmd1 = [ 'cmake "' pathBuildExpressionGraph  '" -G "' compiler ,...
    '" -DTARGET_NAME:STRING="' filename '" -DSDK_DIR:PATH="' SDK_DIR '" -DCPP_DIR:PATH="' cpp_dir ' "'];
system(cmd1)
cmd2 = "cmake --build . --config RelWithDebInfo";
system(cmd2)

%
cd(BIN_DIR)
path_EXE = fullfile(pathBuild, 'RelWithDebInfo', [filename '.exe']);
system(path_EXE)

%% Part 2: build external function (i.e., build .dll).

fooName = 'foo.py';
pathBuildExternalFunction = fullfile(pathMain, 'buildExternalFunction');
path_external_filename_foo = fullfile(BIN_DIR, fooName);
path_external_functions_filename_build = fullfile(pathMain,'build-ExternalFunction',filename);
path_external_functions_filename_install = fullfile(pathMain, 'install-ExternalFunction',filename);
mkdir(path_external_functions_filename_build);
mkdir(path_external_functions_filename_install);
mkdir(pathBuildExternalFunction);
copyfile(path_external_filename_foo, fullfile(pathBuildExternalFunction,fooName));
cd(pathMain)


end