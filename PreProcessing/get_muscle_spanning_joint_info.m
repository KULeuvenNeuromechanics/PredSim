function [muscle_spanning_joint_info] = get_muscle_spanning_joint_info(S,osim_path,model_info)
% --------------------------------------------------------------------------
% get_muscle_spanning_joint_info
%   Find out which muscles span wich joint, thus interacts with its
%   associated coordinates.
% 
% INPUT:
%   - S -
%   * setting structure S
%
%   - osim_path -
%   * path to the OpenSim model file (.osim)
% 
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - muscle_spanning_joint_info -
%   * table with a column for each coordinate and a row for each muscle. 1
%   means this muscle and coordinate interact, 0 means they don't
%
% Original author: Lars D'Hondt
% Original date: 05/April/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------

%%

% names
coordinate_names = model_info.ExtFunIO.coord_names.all;
muscle_names = model_info.muscle_info.muscle_names;

% sizes and indices
n_muscle = model_info.muscle_info.NMuscle;
n_coord = length(coordinate_names);
n_data_points = 5;

% get dummy motion
Qs = generate_dummy_motion(S,model_info,n_data_points);

%% Initialise model
import org.opensim.modeling.*;
model = Model(osim_path);
s = model.initSystem;
% Get state vector
state_vars = model.getStateVariableValues(s);
% Get set of muscles
muscles = model.getMuscles();

%% Find combinations
% Set state vector to 0
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);

dM = zeros(n_data_points,n_muscle,n_coord);

% Loop through dummy states
for j=1:n_data_points
    % Set each coordinate value
    for i=1:n_coord
        state_vars.set((i-1)*2,Qs(j,i));
    end
    model.setStateVariableValues(s,state_vars);
    model.realizePosition(s);

    % Loop over muscles
    for m=1:n_muscle
        muscle_m = muscles.get(muscle_names{m});

        % Get moment arm for each joint
        for i=1:n_coord
            dM(j,m,i) = muscle_m.computeMomentArm(s,model.getCoordinateSet().get(i-1));
        end
    end

end

muscle_spanning_joint_info = squeeze(sum(abs(dM), 1));
muscle_spanning_joint_info(muscle_spanning_joint_info<=0.0001 & muscle_spanning_joint_info>=-0.0001) = 0;
muscle_spanning_joint_info(muscle_spanning_joint_info~=0) = 1;

