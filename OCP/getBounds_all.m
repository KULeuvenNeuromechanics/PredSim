% This script provides bounds and scaling factors for the design variables.
% The bounds on the joint variables are informed by experimental data.
% The bounds on the remaining variables are fixed.
% The bounds are scaled such that the upper/lower bounds cannot be
% larger/smaller than 1/-1.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% Adapted: Lars D'Hondt
% Date: 01 dec 2021
% Adapted: Dhruv Gupta
% Date: 17 jan 2022
%
%--------------------------------------------------------------------------
function [bounds,scaling] = getBounds_all(Qs,model_info,S)

% Get the names of the coordinates
coordinate_names = fieldnames(model_info.ExtFunIO.coordi);
NCoord = length(coordinate_names);
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

%% get IK-besed bounds from spline
% The extreme values are selected as upper/lower bounds, which are then
% further extended.

% prepare index arrays for later use
idx_mtp = [];
idx_shoulder_flex = [];
for i = [model_info.ExtFunIO.jointi.ground_pelvis model_info.ExtFunIO.jointi.back]
    coordinate = coordinate_names{i};
    coord_idx = model_info.ExtFunIO.coordi.(coordinate);
    spline_idx = find(strcmp(Qs.colheaders(1,:),coordinate));
    % Qs
    bounds.Qs.upper(coord_idx) = max(Qs_spline.data(:,spline_idx));
    bounds.Qs.lower(coord_idx) = min(Qs_spline.data(:,spline_idx));
    % Qdots
    bounds.Qdots.upper(coord_idx) = max(Qdots_spline.data(:,spline_idx));
    bounds.Qdots.lower(coord_idx) = min(Qdots_spline.data(:,spline_idx));
    % Qdotdots
    bounds.Qdotdots.upper(coord_idx) = max(Qdotdots_spline.data(:,spline_idx));
    bounds.Qdotdots.lower(coord_idx) = min(Qdotdots_spline.data(:,spline_idx));
end
for i=[model_info.ExtFunIO.jointi.leg_r model_info.ExtFunIO.jointi.arm_r]
    coordinate_r = coordinate_names{i};
    coordinate_l = [coordinate_r(1:end-2) '_l'];
    coord_idx_r = model_info.ExtFunIO.coordi.(coordinate_r);
    coord_idx_l = model_info.ExtFunIO.coordi.(coordinate_l);
    coord_idx = [coord_idx_r coord_idx_l];
    spline_idx_r = find(strcmp(Qs.colheaders(1,:),coordinate_r));
    spline_idx_l = find(strcmp(Qs.colheaders(1,:),coordinate_l));
    spline_idx = [spline_idx_r spline_idx_l];
    % Qs
    bounds.Qs.upper(coord_idx) = max(max(Qs_spline.data(:,spline_idx)));
    bounds.Qs.lower(coord_idx) = min(min(Qs_spline.data(:,spline_idx)));
    % Qdots
    bounds.Qdots.upper(coord_idx) = max(max(Qdots_spline.data(:,spline_idx)));
    bounds.Qdots.lower(coord_idx) = min(min(Qdots_spline.data(:,spline_idx)));
    % Qdotdots
    bounds.Qdotdots.upper(coord_idx) = max(max(Qdotdots_spline.data(:,spline_idx)));
    bounds.Qdotdots.lower(coord_idx) = min(min(Qdotdots_spline.data(:,spline_idx)));
    
    % save indices for later use
    if contains(coordinate_r,'mtp')
        idx_mtp = coord_idx;
    end
    if contains(coordinate_r,'arm_add')
        idx_sh_add = coord_idx;
        rec_uw_sh_add = bounds.Qs.upper(coord_idx);
    elseif contains(coordinate_r,'arm_rot')
        idx_sh_rot = coord_idx;
        rec_uw_sh_rot = bounds.Qs.upper(coord_idx);
    elseif contains(coordinate_r,'elbow_flex')
        idx_elb = coord_idx;
    elseif contains(coordinate_r,'arm_flex')
        idx_shoulder_flex = coord_idx;
    end
end

%% extend IK-based bounds

% The bounds are extended by twice the absolute difference between upper
% and lower bounds.
Qs_range = bounds.Qs.upper - bounds.Qs.lower;
bounds.Qs.lower = bounds.Qs.lower - 2*Qs_range;
bounds.Qs.upper = bounds.Qs.upper + 2*Qs_range;

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
bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(4)) = 2;  
bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(4)) = 0;
% Pelvis_ty
bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(5)) = S.subject.IG_pelvis_y*1.2;
bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(5)) = S.subject.IG_pelvis_y*0.8;
% Pelvis_tz
bounds.Qs.upper(model_info.ExtFunIO.jointi.floating_base(6)) = 0.1;
bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(6)) = -0.1;
% Mtp
bounds.Qs.upper(idx_mtp) = 1.05;
bounds.Qs.lower(idx_mtp) = -0.5;
bounds.Qdots.upper(idx_mtp) = 13;
bounds.Qdots.lower(idx_mtp) = -13;
bounds.Qdotdots.upper(idx_mtp) = 500;
bounds.Qdotdots.lower(idx_mtp) = -500;

bounds.Qs.upper(idx_sh_add) = rec_uw_sh_add;
bounds.Qs.upper(idx_sh_rot) = rec_uw_sh_rot;
bounds.Qs.lower(idx_elb) = 0;

% We adjust some bounds when we increase the speed to allow for the
% generation of running motions.
if S.subject.v_pelvis_x_trgt > 1.33
    % Pelvis tilt
    bounds.Qs.lower(model_info.ExtFunIO.jointi.floating_base(1)) = -20*pi/180;
    % Shoulder flexion
    bounds.Qs.lower(idx_shoulder_flex) = -50*pi/180;
    % Pelvis tx
    bounds.Qdots.upper(model_info.ExtFunIO.jointi.floating_base(4)) = 4;
end

%% Muscle activations
bounds.activation.lower = S.bounds.a.lower*ones(1,NMuscle);
bounds.activation.upper = ones(1,NMuscle);

%% Muscle-tendon forces
bounds.FTtilde.lower = zeros(1,NMuscle);
bounds.FTtilde.upper = 5*ones(1,NMuscle);

%% Time derivative of muscle activations
tact = 0.015;
tdeact = 0.06;
bounds.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tdeact);
bounds.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tact);

%% Time derivative of muscle-tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscle);
bounds.dFTtilde.upper = 1*ones(1,NMuscle);

%% Actuator activations
bounds.activation_actuators.lower = -ones(1,model_info.actuator_info.NActuators);
bounds.activation_actuators.upper = ones(1,model_info.actuator_info.NActuators);

%% Actuator excitations
bounds.excitation_actuators.lower = -ones(1,model_info.actuator_info.NActuators);
bounds.excitation_actuators.upper = ones(1,model_info.actuator_info.NActuators);

%% Final time
bounds.tf.lower = S.bounds.t_final.lower;
bounds.tf.upper = S.bounds.t_final.upper;

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

% Qdotdots
scaling.Qdotdots = max(abs(bounds.Qdotdots.lower),...
    abs(bounds.Qdotdots.upper));
bounds.Qdotdots.lower = (bounds.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0;
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% Arm torque actuators
% Fixed scaling factor
scaling.Actuator_torque = model_info.actuator_info.max_torque;
% Time derivative of muscle activations
% Fixed scaling factor
scaling.vA = 100;
% Muscle activations
scaling.activation = 1;
% Arm activations
scaling.activation_actuators = 1;
% Arm excitations
scaling.excitation_actuators = 1;
% Time derivative of muscle-tendon forces
% Fixed scaling factor
scaling.dFTtilde = 100;
% Muscle-tendon forces
scaling.FTtilde         = max(...
    abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;

%% Hard bounds
% We impose the initial position of pelvis_tx to be 0
bounds.Qs_0.lower = bounds.Qs.lower;
bounds.Qs_0.upper = bounds.Qs.upper;
bounds.Qdots_0.lower = bounds.Qdots.lower;
bounds.Qdots_0.upper = bounds.Qdots.upper;
bounds.Qs_0.lower(model_info.ExtFunIO.jointi.floating_base(4)) = 0;
bounds.Qs_0.upper(model_info.ExtFunIO.jointi.floating_base(4)) = 0;

end