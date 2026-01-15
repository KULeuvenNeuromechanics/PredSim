function [S] = osim2dll(S, osim_path)
% --------------------------------------------------------------------------
% osim2dll
%   This functions uses the OpenSim model to construct a CasADi external
%   function (.dll file) that contains an implicit formulation of the
%   skeletal and contact dynamics. More information can be found here:
%   https://github.com/KULeuvenNeuromechanics/opensimAD
%
% See also generateExternalFunction
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
if ispc
    libextension = '.dll';
elseif isunix
    libextension = '.so';
elseif S.OpenSimADOptions.useSerialisedFunction
    libextension = '.casadi';
end
external_function_dll = fullfile(S.misc.subject_path,['F_', osim_file_name, libextension]);
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

% test 2: do external function files contain all inputs/outputs needed
if extFunOk
    load(external_function_IO,'IO');

    % has field with inputs
    extFunOk = isfield(IO,'input');

    % input forces on bodies
    if extFunOk && ~isempty(S.OpenSimADOptions.input3DBodyForces)
        extFunOk = isfield(IO.input,'Forces');
        for i=1:length(S.OpenSimADOptions.input3DBodyForces)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.input.Forces,S.OpenSimADOptions.input3DBodyForces(i).name);
        end
    end

    % input moments on bodies
    if extFunOk && ~isempty(S.OpenSimADOptions.input3DBodyMoments)
        extFunOk = isfield(IO.input,'Moments');
        for i=1:length(S.OpenSimADOptions.input3DBodyMoments)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.input.Moments,S.OpenSimADOptions.input3DBodyMoments(i).name);
        end
    end

    % output total GRF per foot
    if extFunOk && S.OpenSimADOptions.exportGRFs
        extFunOk = isfield(IO,'GRFs') && isfield(IO.GRFs,'right_total')...
            && isfield(IO.GRFs,'left_total');
    end

    % output GRF per contact sphere
    if extFunOk && S.OpenSimADOptions.exportSeparateGRFs
        if isfield(IO,'GRFs')
            fields_GRF = fieldnames(IO.GRFs);
            extFunOk = any(~strcmp(fields_GRF,'right_total') & ~strcmp(fields_GRF,'left_total'));
        else
            extFunOk = false;
        end
    end

    % output GRM per foot
    if extFunOk && S.OpenSimADOptions.exportGRMs
        extFunOk = isfield(IO,'GRMs') && isfield(IO.GRMs,'right_total')...
            && isfield(IO.GRMs,'left_total');
    end

    % output power per contact sphere
    if extFunOk && S.OpenSimADOptions.exportContactPowers
        extFunOk = isfield(IO,'P_contact_deformation_y');
    end

    % output point positions
    if extFunOk && ~isempty(S.OpenSimADOptions.export3DPositions)
        extFunOk = isfield(IO,'position');
        for i=1:length(S.OpenSimADOptions.export3DPositions)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.position,S.OpenSimADOptions.export3DPositions(i).name);
        end
    end

    % output point velocities
    if extFunOk && ~isempty(S.OpenSimADOptions.export3DVelocities)
        extFunOk = isfield(IO,'velocity');
        for i=1:length(S.OpenSimADOptions.export3DVelocities)
            if ~extFunOk
                break
            end
            extFunOk = isfield(IO.velocity,S.OpenSimADOptions.export3DVelocities(i).name);
        end
    end

    if ~extFunOk
        disp(['   External function does not contain the required inputs and outputs. Removing files:'])
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

if extFunOk && ~S.OpenSimADOptions.always_generate
    disp(['   Using existing external function: '])
    disp(['      ' external_function_dll])
else

    t0 = tic;
    disp('   Creating new external function...')

    if ispc
        addpath(fullfile(S.misc.main_path,'opensimAD','utilities'))
    elseif isunix
        addpath(fullfile(S.misc.main_path,'opensimAD_linux','utilities'))
    end

    outputFilename = ['F_' osim_file_name];

    % Not a setting, to protect users from themselves.
    % Adding expression graphs to calculate 2nd order derivatives for large
    % model (e.g. gait2392) causes size of foo_jac.c file to explode (~10
    % mil lines). Compiler becomes bottleneck.
    secondOrderDerivatives = false;

    try

        if ispc % windows

            generateExternalFunction(osim_path, S.misc.subject_path,...
                S.OpenSimADOptions.jointsOrder, S.OpenSimADOptions.coordinatesOrder,...
                S.OpenSimADOptions.input3DBodyForces, S.OpenSimADOptions.input3DBodyMoments,...
                S.OpenSimADOptions.export3DPositions, S.OpenSimADOptions.export3DOrientations,...
                S.OpenSimADOptions.export3DVelocities, S.OpenSimADOptions.export3DVelocitiesProjGround,...
                S.OpenSimADOptions.exportGRFs, S.OpenSimADOptions.exportGRMs,...
                S.OpenSimADOptions.exportSeparateGRFs, S.OpenSimADOptions.exportContactPowers,...
                outputFilename, S.OpenSimADOptions.compiler,...
                S.OpenSimADOptions.verbose_mode, S.OpenSimADOptions.verify_ID,...
                secondOrderDerivatives, S.OpenSimADOptions.useSerialisedFunction);
        
        else % linux
            generateExternalFunction(osim_path, S.misc.subject_path,...
                S.OpenSimADOptions.jointsOrder, S.OpenSimADOptions.coordinatesOrder,...
                S.OpenSimADOptions.input3DBodyForces, S.OpenSimADOptions.input3DBodyMoments,...
                S.OpenSimADOptions.export3DPositions, S.OpenSimADOptions.export3DVelocities,...
                S.OpenSimADOptions.exportGRFs, S.OpenSimADOptions.exportGRMs,...
                S.OpenSimADOptions.exportSeparateGRFs, S.OpenSimADOptions.exportContactPowers,...
                outputFilename, S.OpenSimADOptions.compiler,...
                S.OpenSimADOptions.verbose_mode, S.OpenSimADOptions.verify_ID,...
                secondOrderDerivatives);
        end

    catch err
        % MATLAB errors in OpenSimAD occur after the actual cause of the
        % error (cmake, compiler, etc.). Print a warning to help the user
        % with troubleshooting the actual problem.
        warn_text = sprintf("OpenSimAD encountered an error.\n");
        if S.OpenSimADOptions.verbose_mode
            warn_text = warn_text + sprintf("Look for error messages in " + ...
                "the text above to identify the cause.\n");
        else
            warn_text = warn_text + ...
                sprintf("Set S.OpenSimADOptions.verbose_mode = true; " + ...
                "and re-run your simulation to get more information.\n");
        end
        warning(warn_text)

        % Still give the MATLAB error incase it has useful info and to stop
        % the simulation
        rethrow(err)
    end

    disp(['   ... external function created (' num2str(toc(t0),'%.2f') ' s)'])

end

S.misc.external_function = ['F_' osim_file_name libextension];

end % end of function
