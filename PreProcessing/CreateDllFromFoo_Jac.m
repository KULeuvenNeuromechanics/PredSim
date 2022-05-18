function [] = CreateDllFromFoo_Jac(filename, cpp_dir, compiler, pathMain, varargin)
%CreateDllFromFoo_Jac Creates a .dll file from foo_jac.c
%   Input arguments:
%       (1) filename: name of the .cpp file
%       (2) cpp_dir: directory .cpp file
%       (3) compiler:  e.g. Visual Studio 15 2017 Win64
%       (4) varargin: 
%           (4.1) verbose_mode: 0 no output, 1 all output

% variable input arguments
verbose_mode = 1;
if ~isempty(varargin)
    verbose_mode = varargin{1};
end

% path information
pathBuildExternalFunction = fullfile(pathMain, 'buildExternalFunction');
path_external_functions_filename_build = fullfile(pathMain,'build-ExternalFunction',filename);
path_external_functions_filename_install = fullfile(pathMain, 'install-ExternalFunction',filename);

% cd to build path
cd(path_external_functions_filename_build);

% run cmake
% compiler = 'Visual Studio 15 2017';
cmd1 = [ 'cmake "' pathBuildExternalFunction  '" -G "' compiler ,...
    '" -DTARGET_NAME:STRING="' filename '" -DINSTALL_DIR:PATH="',...
    path_external_functions_filename_install '"'];
if verbose_mode == 0
    [~,~] = system(cmd1);
else
    system(cmd1);
end
cmd2 = "cmake --build . --config RelWithDebInfo --target install";
if verbose_mode == 0
    [~,~] = system(cmd2);
else
    system(cmd2);
end
cd(pathMain)
copyfile(fullfile(path_external_functions_filename_install, 'bin',[filename '.dll']),...
    fullfile(cpp_dir,[filename '.dll']));

% cleaning up
delete(fullfile(pathBuildExternalFunction, 'foo_jac.c'));
delete(fullfile(pathBuildExternalFunction, 'foo.py'));

% remove specific folders from buildExpressionGraph, build-ExternalFunction
% and install-ExternalFunction (Does not work, I need admin rights for
% this)
% rmdir(path_external_functions_filename_build);
% rmdir(path_external_functions_filename_install);
% pathBuildExpressionGraph = fullfile(pathMain, 'buildExpressionGraph');
% pathBuild = fullfile(pathBuildExpressionGraph,filename);
% rmdir(pathBuild);


end