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
% OUTPUT:
%   - MuscleData -
%   * structure with joint angles and according muscle-tendon lengths and
%       momentarms (input to PolynomialFit.m)
% 
% note: The boundaries of the dummy motion are hard-coded. They should
%       become user inputs in a future update;
%
% Original author: Lars D'Hondt
% Original date: 23/March/2022
%
% Last edit by: 
% Last edit date: 
% --------------------------------------------------------------------------


%% Create dummy motion

% bounds of training dataset (in degrees)
% hard-coded for now, make input later
Bounds(1,:) = {'hip_flex',[-50 50]};
Bounds(2,:) = {'hip_add',[-30 30]};
Bounds(3,:) = {'hip_rot',[-30 30]};
Bounds(4,:) = {'knee',[-90 0]};
Bounds(5,:) = {'ankle',[-30 30]};
Bounds(6,:) = {'subt',[-30 30]};
Bounds(7,:) = {'mtj',[-30 30]};
Bounds(8,:) = {'mtp',[-20 50]};
Bounds(9,:) = {'lumbar_ext',[-30 30]};
Bounds(10,:) = {'lumbar_bend',[-30 30]};
Bounds(11,:) = {'lumbar_rot',[-30 30]};


% names
coordinate_names = model_info.ExtFunIO.coordinate_names;
muscle_names = model_info.muscle_info.muscle_names;

% sizes and indices
n_muscle = model_info.muscle_info.NMuscle;
n_coord = length(coordinate_names);
n_data_points = 5000;
n_data_points = 100; % for debugging only !!!

% default bounds
Q_bounds = [-30;30]*ones(1,n_coord);

% adjust bounds
for j=1:size(Bounds,1)
    idx = find(contains(coordinate_names,Bounds{j,1}));
    if ~isempty(idx)
        Q_bounds(:,idx) = reshape(Bounds{j,2},2,1)*ones(1,length(idx));
    end
end

% construct scale from bounds
Q_scale = diff(Q_bounds);

% generate random joint angles (range 0-1)
Qs = lhsdesign(n_data_points,n_coord);
Qs_test = lhsdesign(5,n_coord);

% scale and offset random joint angles to fit bounds
Qs = Qs.*(ones(n_data_points,1)*Q_scale) + ones(n_data_points,1)*Q_bounds(1,:);
Qs_test = Qs_test.*(ones(5,1)*Q_scale) + ones(5,1)*Q_bounds(1,:);

% base does not have to move for analysis
Qs(:,model_info.ExtFunIO.jointi.floating_base) = 0;
Qs_test(:,model_info.ExtFunIO.jointi.floating_base) = 0;

% angles from degrees to radians
Qs(:,model_info.ExtFunIO.jointi.rotations) = Qs(:,model_info.ExtFunIO.jointi.rotations)*pi/180;
Qs_test(:,model_info.ExtFunIO.jointi.rotations) = Qs_test(:,model_info.ExtFunIO.jointi.rotations)*pi/180;

%% Initialise model
import org.opensim.modeling.*;
model = Model(osim_path);
s = model.initSystem;
% Get state vector
state_vars = model.getStateVariableValues(s);
% Get set of muscles
muscles = model.getMuscles();

%% Find sensitive combinations to evaluate
% Set state vector to 0
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);

dM_test = zeros(5,n_muscle,n_coord);

% Loop through dummy states
for j=1:5
    % Set each coordinate value
    for i=1:n_coord
        state_vars.set((i-1)*2,Qs_test(j,i));
    end
    model.setStateVariableValues(s,state_vars);
    model.realizePosition(s);

    % Loop over muscles
    for m=1:n_muscle
        muscle_m = muscles.get(muscle_names{m});

        % Get moment arm for each joint
        for i=1:n_coord
            dM_test(j,m,i) = muscle_m.computeMomentArm(s,model.getCoordinateSet().get(i-1));
        end
    end

end

muscle_spanning_joint_INFO = squeeze(sum(abs(dM_test), 1));
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO<=0.0001 & muscle_spanning_joint_INFO>=-0.0001) = 0;
muscle_spanning_joint_INFO(muscle_spanning_joint_INFO~=0) = 1;

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
            if muscle_spanning_joint_INFO(m,i)
                dM(j,m,i) = muscle_m.computeMomentArm(s,model.getCoordinateSet().get(i-1));
            end
        end
    end

end


%% Store analysis results
% structure and fieldnames of MuscleData are to be consistent with
% PolynomialFit.m

% coordinate names
MuscleData.dof_names = coordinate_names;

% muscle names
MuscleData.muscle_names = muscle_names;

% joint angles training data
MuscleData.q = Qs;

% muscle-tendon lengths
MuscleData.lMT = lMT;

% moment arms
MuscleData.dM = dM;







