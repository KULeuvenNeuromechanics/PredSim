function [MuscleData] = muscleAnalysisAPI(S,osim_path,model_info)
% --------------------------------------------------------------------------
% muscleAnalysisAPI
%   Analyse the musculoskeletal geometry of the given osim model.
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
%   - run_full_analysis - 
%   * Set to 0 to get muscle_spanning_joint_info as output instead of
%   performing the full analysis. Default value is 1;
%
% OUTPUT:
%   - MuscleData -
%   * structure with joint angles and according muscle-tendon lengths and
%       momentarms (input to PolynomialFit.m)
%
% Original author: Lars D'Hondt
% Original date: 5/April/2022
%
% Last edit by:
% Last edit date:
% --------------------------------------------------------------------------


% names of muscles
muscle_names = model_info.muscle_info.muscle_names;
% number of muscles
n_muscle = model_info.muscle_info.NMuscle;
% number of coordinates
n_coord = length(coordinate_names);
% get senseble muscle-coordinate combinateions to evaluate
muscle_spanning_joint_info = model_info.muscle_info.muscle_spanning_joint_info;

% number of data points
n_data_points = 5000;
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

%% Evaluate muscle-tendon unit lenght and moment arms
% Set state vector to 0
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);

% Initialise matrices for results
lMT = zeros(n_data_points,n_muscle);
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
        % Get MTU length
        lMT(j,m) = muscle_m.getLength(s);
        % Get moment arm for each joint
        for i=1:n_coord
            if muscle_spanning_joint_info(m,i)
                dM(j,m,i) = muscle_m.computeMomentArm(s,model.getCoordinateSet().get(i-1));
            end
        end
    end

end


%% Store analysis results
% structure and fieldnames of MuscleData are to be consistent with PolynomialFit.m

% coordinate names
MuscleData.dof_names = model_info.ExtFunIO.coordinate_names;
% muscle names
MuscleData.muscle_names = muscle_names;
% joint angles training data
MuscleData.q = Qs;
% muscle-tendon lengths
MuscleData.lMT = lMT;
% moment arms
MuscleData.dM = dM;







