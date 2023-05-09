function [S] = osim2dll(S, osim_path)
% --------------------------------------------------------------------------
% osim2dll
%   This functions uses the OpenSim model to construct a CasADi external
%   function (.dll file) that contains an implicit formulation of the
%   skeletal- and contact dynamics. More information can be found here:
%   https://github.com/Lars-DHondt-KUL/opensimAD 
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
% OUTPUT:
%   - S -
%   * setting structure S
% 
% Original author: Lars D'Hondt
% Original date: 9/May/2023
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

% external function files
[~,osim_file_name,~] = fileparts(osim_path);
external_function_dll = fullfile(S.misc.subject_path,['F_', osim_file_name, '.dll']);
external_function_lib = fullfile(S.misc.subject_path,['F_', osim_file_name, '.lib']); % not needed, yet...
external_function_cpp = fullfile(S.misc.subject_path,['F_', osim_file_name, '.cpp']);
external_function_IO  = fullfile(S.misc.subject_path,['F_' osim_file_name '_IO.mat']);

% baseline assumption: there exist suitable external function files
extFunOk = true;

% test 1: do external function files exist
if extFunOk && ~isfile(external_function_dll)
    extFunOk = false;
end
if extFunOk && ~isfile(external_function_IO)
    extFunOk = false;
end

% test 2: existing external function files contain all outputs needed
if extFunOk
    load(external_function_IO,'IO');

    if extFunOk && S.OpenSimADOptions.exportGRFs
        extFunOk = isfield(IO,'GRFs') && isfield(IO.GRFs,'right_foot')...
            && isfield(IO.GRFs,'left_foot');
    end

    if extFunOk && S.OpenSimADOptions.exportSeparateGRFs
        extFunOk = isfield(IO,'GRFs') && isfield(IO.GRFs,'contact_sphere_0');
    end

    if extFunOk && S.OpenSimADOptions.exportGRMs
        extFunOk = isfield(IO,'GRMs') && isfield(IO.GRMs,'right_total')...
            && isfield(IO.GRMs,'left_total');
    end

    if extFunOk && S.OpenSimADOptions.exportContactPowers
        extFunOk = isfield(IO,'P_contact_deformation_y');
    end

    if extFunOk && ~isempty(S.OpenSimADOptions.export3DPositions)
        extFunOk = isfield(IO,'position');
        for i=1:length(S.OpenSimADOptions.export3DPositions)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.position,S.OpenSimADOptions.export3DPositions(i).name);
        end
    end

    if extFunOk && ~isempty(S.OpenSimADOptions.export3DVelocities)
        extFunOk = isfield(IO,'velocity');
        for i=1:length(S.OpenSimADOptions.export3DVelocities)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.velocities,S.OpenSimADOptions.export3DVelocities(i).name);
        end
    end

    if ~extFunOk
        disp(['   External function does not contain the required outputs. Removing files:'])
        disp(['      ' external_function_dll])
        delete(external_function_dll)
        disp(['      ' external_function_IO])
        delete(external_function_IO)
        if isfile(external_function_lib)
            disp(['      ' external_function_lib])
            delete(external_function_lib)
        end
        if isfile(external_function_cpp)
            disp(['      ' external_function_cpp])
            delete(external_function_cpp)
        end
    end
end

if extFunOk
    disp(['   Using existing external function: '])
    disp(['      ' external_function_dll])
else

    t0 = tic;
    disp('   Creating new external function...')

    addpath(fullfile(S.misc.main_path,'opensimAD','utilities'))

    outputFilename = ['F_' osim_file_name];

    generateExternalFunction(osim_path, S.misc.subject_path,...
        S.OpenSimADOptions.jointsOrder, S.OpenSimADOptions.coordinatesOrder,...
        S.OpenSimADOptions.export3DPositions, S.OpenSimADOptions.export3DVelocities,...
        S.OpenSimADOptions.exportGRFs, S.OpenSimADOptions.exportGRMs,...
        S.OpenSimADOptions.exportSeparateGRFs, S.OpenSimADOptions.exportContactPowers,...
        outputFilename, S.OpenSimADOptions.compiler,...
        S.OpenSimADOptions.verbose_mode, S.OpenSimADOptions.verify_ID);

    disp(['   ... external function created (' num2str(toc(t0),'%.2f') ' s)'])

end

S.misc.external_function = ['F_' osim_file_name '.dll'];

end % end of function