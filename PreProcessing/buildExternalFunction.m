function [pathBuildExternalFunction] = buildExternalFunction(filename, cpp_dir,...
    compiler, pathMain,varargin)
% --------------------------------------------------------------------------
% buildExternalFunction
%   Builds external function based on opensim AD modified version of the 
%   OpenSim source code using the installed libraries.
% 
% INPUT:
%   - filename -
%   * name of the cpp file
%
%   - cpp_dir -
%   * folder with cpp file
% 
%   - compiler -
%   * e.g. Visual Studio 15 2017 Win64
%
%   - pathMain -
%   * path to the main directory of the opensimAD folder 
%   (https://github.com/Lars-DHondt-KUL/opensimAD)
%
%   - verbose_mode (optional input) -
%   * 0 no output, 1 all output (default)
% OUTPUT:
%   - model_info -
%   * structure with all the model information based on the OpenSim model
% 
% Original author: Maarten Afschrift
% Original date: 17/May/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


verbose_mode = 1;
if ~isempty(varargin)
    verbose_mode = varargin{1};
end
%% Part 1: create executable from .cpp file

pathBuildExpressionGraph = fullfile(pathMain, 'buildExpressionGraph');
pathBuild = fullfile(pathBuildExpressionGraph,filename);
if ~isfolder(pathBuild)
    mkdir(pathBuild);
end

OpenSimAD_DIR = fullfile(pathMain, 'OpenSimAD-install');
SDK_DIR = fullfile(OpenSimAD_DIR, 'sdk');
BIN_DIR = fullfile(OpenSimAD_DIR, 'bin');

% run cmake to create expression graph
cd(pathBuild)
cmd1 = [ 'cmake "' pathBuildExpressionGraph  '" -G "' compiler ,...
    '" -DTARGET_NAME:STRING="' filename '" -DSDK_DIR:PATH="' SDK_DIR '" -DCPP_DIR:PATH="' cpp_dir ' "'];
if verbose_mode == 0
    [~,~] = system(cmd1);
else
    system(cmd1);
end
cmd2 = "cmake --build . --config RelWithDebInfo";
if verbose_mode == 0
    [~,~] = system(cmd2);
else
    system(cmd2);
end

%
cd(BIN_DIR);
path_EXE = fullfile(pathBuild, 'RelWithDebInfo', [filename '.exe']);
if verbose_mode == 0
    [~,~] =system(path_EXE);
else
    system(path_EXE);
end

%% Part 2: build external function (i.e., build .dll).

fooName = 'foo.py';
pathBuildExternalFunction = fullfile(pathMain, 'buildExternalFunction');
path_external_filename_foo = fullfile(BIN_DIR, fooName);
path_external_functions_filename_build = fullfile(pathMain,'build-ExternalFunction',filename);
path_external_functions_filename_install = fullfile(pathMain, 'install-ExternalFunction',filename);
if isfolder(path_external_functions_filename_build)
    disp(['warning in building dll file: ' path_external_functions_filename_build ...
        ' does already exist. This might cause problems in cmake']);
else
    mkdir(path_external_functions_filename_build);
end
if isfolder(path_external_functions_filename_install)
    disp(['warning in building dll file: ' path_external_functions_filename_install ...
        ' does already exist. This might cause problems in cmake']);
else
    mkdir(path_external_functions_filename_install);
end
if ~isfolder(pathBuildExternalFunction)
    mkdir(pathBuildExternalFunction);
end
copyfile(path_external_filename_foo, fullfile(pathBuildExternalFunction,fooName));
cd(pathMain)


end