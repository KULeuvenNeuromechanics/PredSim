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

for i=1:length(S.orthosis.settings)

    orthosis_settings_i = S.orthosis.settings{i};
    % get Orthosis object
    fun = str2func(orthosis_settings_i.function_name);
    orthosis = fun(init, orthosis_settings_i);

    % run testing methods
    orthosis.setOsimPath(osim_path);
    orthosis.testOsimModel();

    % Add OpenSimAD options required for orthosis
    PointPositions = orthosis.getPointPositions();
    S.OpenSimADOptions.export3DPositions = [S.OpenSimADOptions.export3DPositions(:),...
        PointPositions(:)];
    PointVelocities = orthosis.getPointVelocities();
    S.OpenSimADOptions.export3DVelocities = [S.OpenSimADOptions.export3DVelocities(:),...
        PointVelocities(:)];
    BodyForces = orthosis.getBodyForces();
    S.OpenSimADOptions.input3DBodyForces = [S.OpenSimADOptions.input3DBodyForces(:),...
        BodyForces(:)];
    BodyMoments = orthosis.getBodyMoments();
    S.OpenSimADOptions.input3DBodyMoments = [S.OpenSimADOptions.input3DBodyMoments(:),...
        BodyMoments];
    
    S.orthosis.settings{i}.object = orthosis;

end


end % end of function
