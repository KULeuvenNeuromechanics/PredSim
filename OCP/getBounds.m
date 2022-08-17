function [bounds,scaling] = getBounds(S,model_info)
% --------------------------------------------------------------------------
% getBounds
%   This script provides bounds and scaling factors for the design variables.
%   The bounds on the joint variables are informed by experimental data.
%   The bounds on the remaining variables are fixed.
%   The bounds are scaled such that the upper/lower bounds cannot be
%   larger/smaller than 1/-1.
%   
% INPUT:
%   - S -
%   * setting structure S
%
%   - model_info -
%   * structure with all the model information based on the OpenSim model
%
% OUTPUT:
%   - bounds -
%   * boundaries for all optimisation variables
%
%   - scaling -
%   * scale factors for all optimisation variables
% 
% Original author: Antoine Falisse
% Original date: 12/19/2018
%
%   Edit: When simulating a half gait cycle, the upper bound on pelvis_tx
%   position should be halved, not the bound on its velocity.
% Last edit by: Lars D'Hondt
% Last edit date: 17 aug 2022
% --------------------------------------------------------------------------

% Kinematics file for bounds -- input arguments
% IKfile_bounds = fullfile(S.subject.folder_name, S.subject.IKfile_bounds);
Qs = getIK(S.subject.IK_Bounds,model_info);

% Get the names of the coordinates
coordinate_names = model_info.ExtFunIO.coord_names.all;
NCoord = model_info.ExtFunIO.jointi.nq.all;
NMuscle = model_info.muscle_info.NMuscle;

%% Spline approximation of Qs to get Qdots and Qdotdots
Qs_spline.data = zeros(size(Qs.allfilt));
Qs_spline.data(:,1) = Qs.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allfilt));
Qdots_spline.data(:,1) = Qs.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allfilt));
Qdotdots_spline.data(:,1) = Qs.allfilt(:,1);
for i = 2:size(Qs.allfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1),Qs.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allfilt(:,1),1);
end

%% get IK-based bounds from spline
% The extreme values are selected as upper/lower bounds, which are then
% further extended.

% prepare index arrays for later use
idx_mtp = [];
idx_arms = [model_info.ExtFunIO.jointi.arm_r,model_info.ExtFunIO.jointi.arm_l];
idx_shoulder_flex = [];
idx_elbow = [];
for i = 1:NCoord
    coordinate = coordinate_names{i};
    coord_idx = model_info.ExtFunIO.coordi.(coordinate);
    spline_idx = strcmp(Qs.colheaders(1,:),coordinate);
    % Qs
    bounds.Qs.upper(coord_idx) = max((Qs_spline.data(:,spline_idx)));
    bounds.Qs.lower(coord_idx) = min((Qs_spline.data(:,spline_idx)));
    % Qdots
    bounds.Qdots.upper(coord_idx) = max((Qdots_spline.data(:,spline_idx)));
    bounds.Qdots.lower(coord_idx) = min((Qdots_spline.data(:,spline_idx)));
    % Qdotdots
    bounds.Qdotdots.upper(coord_idx) = max((Qdotdots_spline.data(:,spline_idx)));
    bounds.Qdotdots.lower(coord_idx) = min((Qdotdots_spline.data(:,spline_idx)));

    % save indices for later use
    if contains(coordinate,'mtp')
        idx_mtp(end+1) = coord_idx;
    end
    if find(idx_arms(:)==coord_idx)
        if contains(coordinate,'elbow')
            idx_elbow(end+1) = coord_idx;
        elseif contains(coordinate,'flex')
            idx_shoulder_flex(end+1) = coord_idx;
        end
    end
end

%% extend IK-based bounds
idx_extend = [model_info.ExtFunIO.jointi.floating_base,...
    model_info.ExtFunIO.jointi.leg_r,...
    model_info.ExtFunIO.jointi.leg_l,...
    model_info.ExtFunIO.jointi.torso,...
    idx_elbow,idx_shoulder_flex];

% The bounds are extended by twice the absolute difference between upper
% and lower bounds.
Qs_range = abs(bounds.Qs.upper - bounds.Qs.lower);
bounds.Qs.lower(idx_extend) = bounds.Qs.lower(idx_extend) - 2*Qs_range(idx_extend);
bounds.Qs.upper(idx_extend) = bounds.Qs.upper(idx_extend) + 2*Qs_range(idx_extend);

% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdots_range = abs(bounds.Qdots.upper - bounds.Qdots.lower);
bounds.Qdots.lower = bounds.Qdots.lower - 3*Qdots_range;
bounds.Qdots.upper = bounds.Qdots.upper + 3*Qdots_range;

% The bounds are extended by 3 times the absolute difference between upper
% and lower bounds.
Qdotdots_range = abs(bounds.Qdotdots.upper - bounds.Qdotdots.lower);
bounds.Qdotdots.lower = bounds.Qdotdots.lower - 3*Qdotdots_range;
bounds.Qdotdots.upper = bounds.Qdotdots.upper + 3*Qdotdots_range;

%% manual adjustment
% For several joints, we manually adjust the bounds
% floating base tx
bounds.Qs.upper(model_info.ExtFunIO.jointi.base_forward) = 4;  
bounds.Qs.lower(model_info.ExtFunIO.jointi.base_forward) = 0;
% Pelvis_ty
if length(model_info.ExtFunIO.jointi.floating_base)==6
    bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(5)) = model_info.IG_pelvis_y*1.2;
    bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(5)) = model_info.IG_pelvis_y*0.5;
    
elseif length(model_info.ExtFunIO.jointi.floating_base)==3
    bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(3)) = model_info.IG_pelvis_y*1.2;
    bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(3)) = model_info.IG_pelvis_y*0.5;

end
% Pelvis_tz
bounds.Qs.upper(model_info.ExtFunIO.jointi.base_lateral) = 0.1;
bounds.Qs.lower(model_info.ExtFunIO.jointi.base_lateral) = -0.1;

% Elbow
bounds.Qs.lower(idx_elbow) = 0;
% Mtp
bounds.Qs.upper(idx_mtp) = 1.05;
bounds.Qs.lower(idx_mtp) = -0.5;
bounds.Qdots.upper(idx_mtp) = 13;
bounds.Qdots.lower(idx_mtp) = -13;
bounds.Qdotdots.upper(idx_mtp) = 500;
bounds.Qdotdots.lower(idx_mtp) = -500;

% We adjust some bounds when we increase the speed to allow for the
% generation of running motions.
if S.subject.v_pelvis_x_trgt > 1.33
    % Pelvis tilt
    bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(1)) = -20*pi/180;
    % Shoulder flexion
    bounds.Qs.lower(idx_shoulder_flex) = -50*pi/180;
    % Pelvis tx
    bounds.Qdots.upper(model_info.ExtFunIO.jointi.base_forward) = 8;
end

if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    bounds.Qs.upper(model_info.ExtFunIO.jointi.base_forward) = ...
        bounds.Qs.upper(model_info.ExtFunIO.jointi.base_forward)/2;
end

%% Adjust bounds based on settings
if ~isempty(S.bounds.coordinates)
    [new_lb,new_ub] = unpack_name_value_combinations(S.bounds.coordinates,coordinate_names,[1,1]);
    
    new_lb(model_info.ExtFunIO.jointi.rotations) = new_lb(model_info.ExtFunIO.jointi.rotations)*pi/180;
    new_ub(model_info.ExtFunIO.jointi.rotations) = new_ub(model_info.ExtFunIO.jointi.rotations)*pi/180;
    
    for i = 1:NCoord
        coordinate = coordinate_names{i};
        coord_idx = model_info.ExtFunIO.coordi.(coordinate);
    
        if ~isnan(new_lb(i))
            bounds.Qs.lower(coord_idx) = new_lb(i);
        end
    
        if ~isnan(new_ub(i))
            bounds.Qs.upper(coord_idx) = new_ub(i);
        end
    end

end

%% Adjust bounds to be symmetric
bounds.Qs.upper(model_info.ExtFunIO.symQs.QsInvA) =...
    max(bounds.Qs.upper(model_info.ExtFunIO.symQs.QsInvA),bounds.Qs.upper(model_info.ExtFunIO.symQs.QsInvB));
bounds.Qs.lower(model_info.ExtFunIO.symQs.QsInvA) =...
    min(bounds.Qs.lower(model_info.ExtFunIO.symQs.QsInvA),bounds.Qs.lower(model_info.ExtFunIO.symQs.QsInvB));

bounds.Qdots.upper(model_info.ExtFunIO.symQs.QdotsInvA) =...
    max(bounds.Qdots.upper(model_info.ExtFunIO.symQs.QdotsInvA),bounds.Qdots.upper(model_info.ExtFunIO.symQs.QdotsInvB));
bounds.Qdots.lower(model_info.ExtFunIO.symQs.QdotsInvA) =...
    min(bounds.Qdots.lower(model_info.ExtFunIO.symQs.QdotsInvA),bounds.Qdots.lower(model_info.ExtFunIO.symQs.QdotsInvB));

bounds.Qdotdots.upper(model_info.ExtFunIO.symQs.QdotsInvA) =...
    max(bounds.Qdotdots.upper(model_info.ExtFunIO.symQs.QdotsInvA),bounds.Qdotdots.upper(model_info.ExtFunIO.symQs.QdotsInvB));
bounds.Qdotdots.lower(model_info.ExtFunIO.symQs.QdotsInvA) =...
    min(bounds.Qdotdots.lower(model_info.ExtFunIO.symQs.QdotsInvA),bounds.Qdotdots.lower(model_info.ExtFunIO.symQs.QdotsInvB));



%% Muscle activations
bounds.a.lower = S.bounds.a.lower*ones(1,NMuscle);
bounds.a.upper = ones(1,NMuscle);

%% Muscle-tendon forces
bounds.FTtilde.lower = zeros(1,NMuscle);
bounds.FTtilde.upper = 5*ones(1,NMuscle);

%% Time derivative of muscle activations
bounds.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tdeact);
bounds.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*model_info.muscle_info.tact);

%% Time derivative of muscle-tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscle);
bounds.dFTtilde.upper = 1*ones(1,NMuscle);

%% Torque actuator activations
bounds.a_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds.a_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);

%% Torque actuator excitations
bounds.e_a.lower = -ones(1,model_info.ExtFunIO.jointi.nq.torqAct);
bounds.e_a.upper = ones(1,model_info.ExtFunIO.jointi.nq.torqAct);

%% Final time
bounds.tf.lower = S.bounds.t_final.lower;
if strcmp(S.misc.gaitmotion_type,'HalfGaitCycle')
    bounds.tf.upper = S.bounds.t_final.upper/2;
else
    bounds.tf.upper = S.bounds.t_final.upper;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Scaling
% Qs
scaling.Qs      = max(abs(bounds.Qs.lower),abs(bounds.Qs.upper));
bounds.Qs.lower = (bounds.Qs.lower)./scaling.Qs;
bounds.Qs.upper = (bounds.Qs.upper)./scaling.Qs;
% Qdots
scaling.Qdots      = max(abs(bounds.Qdots.lower),abs(bounds.Qdots.upper));
bounds.Qdots.lower = (bounds.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds.Qdots.upper)./scaling.Qdots;
% Qs and Qdots are intertwined
bounds.QsQdots.lower = zeros(1,2*NCoord);
bounds.QsQdots.upper = zeros(1,2*NCoord);
bounds.QsQdots.lower(1,1:2:end) = bounds.Qs.lower;
bounds.QsQdots.upper(1,1:2:end) = bounds.Qs.upper;
bounds.QsQdots.lower(1,2:2:end) = bounds.Qdots.lower;
bounds.QsQdots.upper(1,2:2:end) = bounds.Qdots.upper;
scaling.QsQdots                 = zeros(1,2*NCoord);
scaling.QsQdots(1,1:2:end)      = scaling.Qs ;
scaling.QsQdots(1,2:2:end)      = scaling.Qdots ;
% Qdotdots
scaling.Qdotdots = max(abs(bounds.Qdotdots.lower),abs(bounds.Qdotdots.upper));
bounds.Qdotdots.lower = (bounds.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0;
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% Torque actuator torque actuators
scaling.ActuatorTorque = struct_array_to_double_array(model_info.actuator_info.parameters,'max_torque');
% Time derivative of muscle activations
% Fixed scaling factor
scaling.vA = 100;
% Muscle activations
scaling.a = 1;
% Torque actuator activations
scaling.a_a = 1;
% Torque actuator excitations
scaling.e_a = 1;
% Time derivative of muscle-tendon forces
% Fixed scaling factor
scaling.dFTtilde = 100;
% Muscle-tendon forces
scaling.FTtilde         = max(abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;

%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds.Qs_0.lower = bounds.Qs.lower;
bounds.Qs_0.upper = bounds.Qs.upper;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_forward) = 0;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.base_lateral) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.base_lateral) = 0;

end
