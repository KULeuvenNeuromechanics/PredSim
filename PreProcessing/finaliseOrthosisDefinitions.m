function [S] = finaliseOrthosisDefinitions(S, osim_path)
% --------------------------------------------------------------------------
% finaliseOrthosisDefinitions
%   Wraps the user-defined orthosis functions. 
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%
% OUTPUT:
%   - S -
%   * setting structure S
%
% 
% Original author: Lars D'Hondt
% Original date: 5/January/2024
% --------------------------------------------------------------------------

% create struct to initialise orthosis
init.Nmesh = S.solver.N_meshes;
init.osimPath = osim_path;

statecounter = uint32(0);
controlcounter = uint32(0);
% first loop: assemble cellstring of all orthosis states and controls
for i=1:length(S.orthosis.settings)
    orthosis_settings_i = S.orthosis.settings{i};

    % get Orthosis object
    fun = str2func(orthosis_settings_i.function_name);
    orthosis = fun(init, orthosis_settings_i);    
    S.orthosis.settings{i}.object = orthosis;

    [~,~,~,~,~,~,meta_arg,~] = getArgRes(orthosis);

    isOptivar = strcmp({meta_arg.type}, 'optivar');

    isX = strcmp({meta_arg.subtype}, 'x') & isOptivar;
    if any(isX)
        S.orthosis.settings{i}.states.names = {meta_arg(isX).name};
        S.orthosis.settings{i}.states.bounds_nsc = {meta_arg(isX).bounds_nsc};
        S.orthosis.settings{i}.states.bounds = {meta_arg(isX).bounds};
        statecounter = statecounter + uint32(1);
    else
        S.orthosis.settings{i}.states.names = {};
        S.orthosis.settings{i}.states.bounds_nsc = {};
        S.orthosis.settings{i}.states.bounds = {};
    end

    isU = strcmp({meta_arg.subtype}, 'u') & isOptivar;
    if any(isU)
        S.orthosis.settings{i}.controls.names = {meta_arg(isU).name};
        S.orthosis.settings{i}.controls.bounds_nsc = {meta_arg(isU).bounds_nsc};
        S.orthosis.settings{i}.controls.bounds = {meta_arg(isU).bounds};
        controlcounter = controlcounter + uint32(1);
    else
        S.orthosis.settings{i}.controls.names = {};
        S.orthosis.settings{i}.controls.bounds_nsc = {};
        S.orthosis.settings{i}.controls.bounds = {};
    end
         

end

% Convert cell array of structs -> struct array
settings = [S.orthosis.settings{:}];

% Concatenate controlNames across all structs
if ~isempty(settings)
    temp = [settings.controls];
    S.orthosis.controlNames_all = [temp.names]; 
    S.orthosis.controlBounds_nsc_all = [temp.bounds_nsc];
    S.orthosis.controlBounds_all = [temp.bounds];
    S.orthosis.Ncontrols_all = controlcounter;
    
    % Concatenate stateNames across all structs
    temp = [settings.states];
    %S.orthosis.Nstates_all = sum()
    S.orthosis.stateNames_all = [temp.names];
    S.orthosis.stateBounds_nsc_all = [temp.bounds_nsc];
    S.orthosis.stateBounds_all = [temp.bounds];
    S.orthosis.Nstates_all = statecounter;

    for i=1:length(S.orthosis.settings)
    
        orthosis = S.orthosis.settings{i}.object;
    
        orthosis.createCasadiFunction();
        % run testing methods
    %     orthosis.setOsimPath(osim_path);
    %     orthosis.testOsimModel();
    
        % Add OpenSimAD options required for orthosis
        PointPositions = orthosis.getPointPositions();
        S.OpenSimADOptions.export3DPositions = [S.OpenSimADOptions.export3DPositions(:);...
            PointPositions(:)];
        PointVelocities = orthosis.getPointVelocities();
        S.OpenSimADOptions.export3DVelocities = [S.OpenSimADOptions.export3DVelocities(:);...
            PointVelocities(:)];
        BodyForces = orthosis.getBodyForces();
        S.OpenSimADOptions.input3DBodyForces = [S.OpenSimADOptions.input3DBodyForces(:);...
            BodyForces(:)];
        BodyMoments = orthosis.getBodyMoments();
        S.OpenSimADOptions.input3DBodyMoments = [S.OpenSimADOptions.input3DBodyMoments(:);...
            BodyMoments(:)];
    
        S.orthosis.settings{i}.object = orthosis;
        
    end
end

%% Remove duplicate entries
fields = ["export3DPositions","export3DVelocities","input3DBodyForces","input3DBodyMoments"];
for j=fields
    
    segments = S.OpenSimADOptions.(j);
    
    if isempty(segments)
        continue
    end
    
    isUnique = true(size(segments));

    for ii = 1:length(segments)-1
        for jj = ii+1:length(segments)
            if isequal(segments(ii),segments(jj))
                isUnique(jj) = false;
                break;
            end
        end
    end
    
    segments(~isUnique) = [];

    S.OpenSimADOptions.(j) = segments;
end

end % end of function
