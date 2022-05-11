function [S] = osim2dll(S,osim_path)
%
% Matlab shell for the python workflow
%
% 
% Author: Lars D'Hondt
%
% Date: 07/March/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to the opensim-AD submodule
pathFunction = mfilename('fullpath');
[pathPreProcessing, ~, ~] = fileparts(pathFunction);
[pathMain, ~] = fileparts(pathPreProcessing);
PathOsimAD = fullfile(pathMain,'opensimAD');

% extract name of opensim model file
[~,osim_file_name,~] = fileparts(osim_path);

% external function name
S.misc.external_function = ['F_' osim_file_name '.dll'];

if ~isfile(fullfile(S.misc.subject_path,S.misc.external_function))

    % list with all input arguments to create the .dll file
    pathOpenSimModel = osim_path;
    outputFilename = ['F_' osim_file_name];
    compiler = S.Cpp2Dll.compiler;
    export3DSegmentOrigins = {'calcn_r', 'calcn_l', 'femur_r', 'femur_l', 'hand_r',...
                          'hand_l', 'tibia_r', 'tibia_l', 'toes_r', 'toes_l'};
    jointsOrder = [];
    coordinatesOrder = [];
    exportGRFs = true;
    exportSeparateGRFs = true;
    exportGRMs = true;
    exportContactPowers = false;
    CppDir = S.misc.subject_path;
    CppName = S.misc.external_function;

    % Create the .cpp file using an executable
    if isempty(S.Cpp2Dll.PathCpp2Dll_Exe)
        error(['You should specify the path to the executables to convert', ...
            'osim file to an .dll file. You can download the executable here: ', ...
            'https://drive.google.com/file/d/1FqlQR1NkUfbRby8K1Kj9e65rF-q_GDkm/view?usp=sharing']);
    PathCpp2Dll_Exe = S.Cpp2Dll.PathCpp2Dll_Exe;

    % create folders if needed
    if ~isfolder(CppDir)
        mkdir(CppDir);
    end

    % create opensim bin folder if needed (quick fix for issues with use of
    % python api in exectutable.
%     if ~isfolder('C:\OpenSim 4.3\bin')
%         mkdir('C:\OpenSim 4.3\bin');
%         copyfile(fullfile(pwd,'OsimBinFolder','bin'), 'C:\OpenSim 4.3\bin');
%     end

    % create a mat file with all input information
    MatFileExoInfo = fullfile(pwd,'MatFileExoInfo.mat');
    pathOut = CppDir;
    save(MatFileExoInfo,'pathOpenSimModel','outputFilename','export3DSegmentOrigins',...
        'jointsOrder','coordinatesOrder','exportGRFs','exportSeparateGRFs',...
        'exportGRMs','exportContactPowers','pathOut');

    % run the exectuable
    MainPath = pwd;
    cd(PathCpp2Dll_Exe);
    command = ['osimtocppexe.exe "' MatFileExoInfo '"'];
    disp('Converting .osim file to .cpp file to solve ID');
    system(command);

    % load .mat file with input and outputs
    load(fullfile(CppDir,[CppName '_IO.mat']),'IO');

    % Create foo file
    cd(MainPath);
    disp('Using cmake to create project for .cpp file');
    [fooPath] = buildExternalFunction(CppName, CppDir, compiler, PathOsimAD);

    % create foo_jac (temporary solution using executable);
    disp('convert foo.py to foo_jac.c using casadi');
    cd(PathCpp2Dll_Exe)
    command = ['GenF.exe ' fooPath ' ' num2str(IO.nCoordinates*3)]; % ToDo: Add nInputes to IO structure
    system(command);
    cd(MainPath);

    % createDll
    disp('Create .dll file using cmake');
    CreateDllFromFoo_Jac(CppName, CppDir, compiler, PathOsimAD);
    cd(MainPath);

    % delete the matfile with all input information again
    delete(MatFileExoInfo);
end


