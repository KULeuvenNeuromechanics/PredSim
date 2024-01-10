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
